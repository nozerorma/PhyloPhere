# perm_stream.py — Replay of the permuted labelings over one gene, in process, scored one base cycle at a time.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
Streamed replay of the permutation null for one gene.

The CAAS kernel of `ct perm-replay` (VectorizedPermReplay, subworkflows/CT/local/modules/perm_replay.py) is called in
the worker itself. The alignment positions are encoded once; the labelings are replayed in ranges of whole base cycles,
and each base cycle is turned into discovery entries, scored against the gene's ASR context, pooled over its
hypotheses and reduced to per-position records before the next one is built. What the worker holds therefore depends
on the labelings of a range and on the records it keeps, not on the discovery rows of the whole gene, and no
discovery export is written.

The entries are built from the same rows the export writes (perm_replay.perm_discovery_rows) and the cycles come out in
the same order as a replay of the export (cycle tags in text order, the hypotheses of a cycle in text order), so the
detail rows equal those of the export route.

Imported by: src/utils/gene_wrapper.py (process_all_genes_perms)
Inputs: a labelings file, a gene's alignment and its ASR context, the replay thresholds of `ct perm-replay`
Outputs: per gene, the pooled per-position records of every cycle with a detection (all sides, and the strict pooling)
"""

import argparse
import logging
import shlex
import sys
import time
from dataclasses import dataclass, fields
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Tuple

from src.core.driver import pool_labelings, score_labelings
from src.core.labelings import base_cycle, read_pss
from src.data.models import CAASPosition
from src.utils.amino import normalize_amino_list

logger = logging.getLogger(__name__)


# ── Replay thresholds ─────────────────────────────────────────────────────────

# Option names, defaults and meaning are those of `ct perm-replay` (subworkflows/CT/local/ct), whose arguments the
# Nextflow task writes to one file; the stream reads that file so the two routes take the same thresholds.
_THRESHOLD_OPTIONS = [
    ("--patterns", "patterns", "1,2,3"),
    ("--max_bg_gaps", "max_bg_gaps", "NO"),
    ("--max_fg_gaps", "max_fg_gaps", "NO"),
    ("--max_gaps", "max_gaps", "NO"),
    ("--max_gaps_per_position", "max_gaps_per_position", "0.5"),
    ("--max_bg_miss", "max_bg_miss", "NO"),
    ("--max_fg_miss", "max_fg_miss", "NO"),
    ("--max_miss", "max_miss", "NO"),
    ("--max_conserved", "max_conserved", "0"),
    ("--min_non_overlapping", "min_non_overlapping", "2"),
]


@dataclass(frozen=True)
class ReplaySettings:
    """Where the CAAS kernel lives, the design the alignment is sliced against and the replay thresholds."""
    ct_dir: str                          # directory holding the `modules` package of `ct`
    design: str                          # trait file, or directory of traitfile_H*.tab, of the run (alignment slicing)
    ali_format: str = "clustal"
    patterns: str = "1,2,3"
    max_bg_gaps: str = "NO"
    max_fg_gaps: str = "NO"
    max_gaps: str = "NO"
    max_gaps_per_position: str = "0.5"
    max_bg_miss: str = "NO"
    max_fg_miss: str = "NO"
    max_miss: str = "NO"
    max_conserved: str = "0"
    min_non_overlapping: str = "2"
    miss_pair: bool = False
    chunk_labelings: Optional[int] = None   # labelings per replay range; None = the kernel's memory budget

    @classmethod
    def from_args(cls, tokens: List[str], ct_dir: str, design: str, ali_format: str,
                  chunk_labelings: Optional[int] = None) -> "ReplaySettings":
        """Settings from `ct perm-replay` style tokens (`--patterns 1,2,3 --miss_pair --max_conserved 2 ...`)."""
        parser = argparse.ArgumentParser(add_help=False)
        for flag, dest, default in _THRESHOLD_OPTIONS:
            parser.add_argument(flag, dest=dest, default=default)
        parser.add_argument("--miss_pair", dest="miss_pair", action="store_true")
        ns = parser.parse_args(tokens)
        values = {f.name: getattr(ns, f.name) for f in fields(cls) if hasattr(ns, f.name)}
        return cls(ct_dir=str(ct_dir), design=str(design), ali_format=ali_format, chunk_labelings=chunk_labelings, **values)

    @classmethod
    def from_args_file(cls, path, ct_dir: str, design: str, ali_format: str,
                       chunk_labelings: Optional[int] = None) -> "ReplaySettings":
        return cls.from_args(shlex.split(Path(path).read_text()), ct_dir, design, ali_format, chunk_labelings)


# ── The labelings and the CAAS kernel of one block ───────────────────────────


def _labelings_file(resample_dir: str) -> Path:
    """The one labelings file of a block directory (`fop_labelings.tab`, else `resample_*.tab`), or the file itself."""
    p = Path(resample_dir)
    if p.is_file():
        return p
    tabs = sorted(p.glob("fop_labelings.tab")) or sorted(p.glob("resample_*.tab"))
    if len(tabs) != 1:
        raise ValueError(f"expected one labelings file under {p}, found {[t.name for t in tabs]}")
    return tabs[0]


class ReplayContext:
    """The labelings of a block, parsed once, and the kernel modules: shared by every gene a worker replays."""

    def __init__(self, settings: ReplaySettings, resample_dir: str, cycles: Optional[List[str]] = None,
                 fop_pairs_file: Optional[str] = None):
        ct_dir = str(Path(settings.ct_dir).resolve())
        if ct_dir not in sys.path:
            sys.path.insert(0, ct_dir)
        from modules import perm_replay as pr
        from modules.perm_replay_io import simtrait_revive
        from modules.runslice import runslice

        self.settings = settings
        self.pr = pr
        self._runslice = runslice
        cfg = simtrait_revive(str(_labelings_file(resample_dir)))
        if cycles is not None:
            wanted = set(cycles)
            cfg.alltraits = [t for t in cfg.alltraits if t in wanted]
        self.cfg = cfg
        # (foreground, background) of every labeling, in the order of the file: what the scoring takes its design from
        self.labelings: Dict[str, Tuple[List[str], List[str]]] = {
            t: (cfg.trait2fg[t], cfg.trait2bg[t]) for t in cfg.alltraits if t in cfg.trait2fg and t in cfg.trait2bg}
        self.fop_pairs = read_pss(fop_pairs_file, labelings=self.labelings) if fop_pairs_file else None

    def slice_alignment(self, alignment_path: str):
        """The alignment as `ct perm-replay` slices it: the columns with at most the allowed share of gaps."""
        s = self.settings
        options = SimpleNamespace(
            single_alignment=str(alignment_path), ali_format=s.ali_format, config_file=s.design,
            max_fg_gaps_string=s.max_fg_gaps, max_bg_gaps_string=s.max_bg_gaps,
            max_fg_miss_string=s.max_fg_miss, max_bg_miss_string=s.max_bg_miss,
            max_gaps_pos_string=s.max_gaps_per_position)
        return self._runslice(options)


# ── Discovery entries and per-position records ───────────────────────────────


def entry_from_fields(fields_: List[str]) -> CAASPosition:
    """The CAASPosition of one perm_discovery row (cycle, gene, position, caas, pattern[, is_conserved_meta, conserved_pair])."""
    pos0 = int(fields_[2].strip())
    caas = fields_[3]
    parts_caas = caas.split("/") if caas else []
    trait1 = normalize_amino_list(list(parts_caas[0])) if len(parts_caas) == 2 else []
    trait0 = normalize_amino_list(list(parts_caas[1])) if len(parts_caas) == 2 else []
    meta = len(fields_) > 6
    cp = fields_[6].strip() if meta else ""
    if cp and ":" in cp.split(",")[0]:
        cp = cp.split(":", 1)[1]
    return CAASPosition(
        position=pos0,
        position_one_based=pos0 + 1,
        tag=f"POS{pos0}",
        caas=caas,
        trait1_aa=trait1,
        trait0_aa=trait0,
        is_conserved_meta=bool(set(trait1) & set(trait0)),
        conserved_pair=cp,
    )


class SlimRecord(NamedTuple):
    """What the detail rows read from a pooled per-position record."""
    position: int
    asr_path_score: float
    side: str
    n_hypotheses: Any
    n_species: Any


def _slim(records: List[Any]) -> List[SlimRecord]:
    return [SlimRecord(r.position, r.asr_path_score, r.side, r.n_hypotheses, r.n_species) for r in records]


def _keep_records(cycle: str, records: List[Any]) -> List[SlimRecord]:
    """The default reduction of a cycle: its records, slimmed."""
    return _slim(records)


# ── One gene ─────────────────────────────────────────────────────────────────


def stream_gene(
    gene: str,
    alignment_path: str,
    replay: ReplayContext,
    load_ctx: Callable[[], Optional[Dict[str, Any]]],
    posterior_threshold: float,
    strict: bool = False,
    reduce: Callable[[str, List[Any]], Any] = _keep_records,
    reduce_strict: Optional[Callable[[str, List[Any]], Any]] = None,
) -> Optional[Tuple[List[Tuple[str, Any]], List[Tuple[str, Any]]]]:
    """Replay every labeling of `replay` over one alignment and score the cycles that have a detection.

    Returns (pooled, strict_pooled): [(cycle or base cycle, reduce(cycle, per-position records))] in cycle-tag order,
    the second list empty unless `strict`; the strict cycles go through `reduce_strict` (default: `reduce`). `reduce` turns the pooled records of a cycle into what the caller keeps (the
    slimmed records by default); a cycle for which it returns None is dropped. ([], []) when no labeling detects
    anything, and None when the gene's ASR context could not be loaded (`load_ctx` is called once, at the first cycle
    with a detection, so a gene without any does not pay for it).
    """
    pr, s = replay.pr, replay.settings
    t_start = time.perf_counter()
    sliced = replay.slice_alignment(alignment_path)
    positions = list(sliced.d)
    vb = pr.VectorizedPermReplay(replay.cfg, sliced.species)
    posdicts = pr._hits_by_position(positions, sliced.genename)
    max_conserved, min_non = int(s.max_conserved), int(s.min_non_overlapping)
    t_setup = time.perf_counter() - t_start

    ctx: Optional[Dict[str, Any]] = None
    t_kernel = t_entries = t_score = t_ctx = 0.0
    n_hits = n_cycles = 0
    pooled: List[Tuple[str, Any]] = []
    pooled_strict: List[Tuple[str, Any]] = []

    chunks = vb.count_chunks(
        positions, sliced.genename, s.max_fg_gaps, s.max_bg_gaps, s.max_gaps, s.max_fg_miss, s.max_bg_miss, s.max_miss,
        max_conserved, s.patterns, b_chunk=s.chunk_labelings, whole_base_cycles=True, collect_hits=True,
        miss_pair=s.miss_pair, min_non_overlapping=min_non)
    while True:
        t0 = time.perf_counter()
        try:
            _start, _end, _counts, hits = next(chunks)
        except StopIteration:
            break
        t_kernel += time.perf_counter() - t0

        # the hits of a range, regrouped by base cycle (position order is kept inside each)
        by_base: Dict[str, Dict[str, List[str]]] = {}
        for position_name, traits in hits.items():
            for trait in traits:
                by_base.setdefault(base_cycle(trait), {}).setdefault(position_name, []).append(trait)

        for base, base_hits in by_base.items():
            t0 = time.perf_counter()
            entries_by_tag: Dict[str, List[CAASPosition]] = {}
            disagreements = []
            for trait, posnum, row, agrees in pr.perm_discovery_rows(
                    replay.cfg, sliced.genename, positions, base_hits, max_conserved, min_non, posdicts):
                if not agrees:
                    disagreements.append((trait, sliced.genename, posnum))
                entries_by_tag.setdefault(trait, []).append(entry_from_fields(row))
                n_hits += 1
            if disagreements:
                raise pr.KernelDisagreement(
                    f"vectorized kernel and check_pattern disagree on {len(disagreements)} hit(s) "
                    f"(labeling, gene, position), first: {disagreements[:5]}")
            t_entries += time.perf_counter() - t0

            if ctx is None:
                t0 = time.perf_counter()
                ctx = load_ctx()
                t_ctx = time.perf_counter() - t0
                if ctx is None:
                    return None

            t0 = time.perf_counter()
            results = score_labelings(ctx, gene, sorted(entries_by_tag), replay.labelings, entries_by_tag, posterior_threshold)
            if results:
                n_cycles += 1
                for cyc, recs in pool_labelings(results, replay.fop_pairs, replay.labelings):
                    kept = reduce(cyc, _slim(recs))
                    if kept is not None:
                        pooled.append((cyc, kept))
                if strict:
                    for cyc, recs in pool_labelings(results, replay.fop_pairs, replay.labelings, strict=True):
                        kept = (reduce_strict or reduce)(cyc, _slim(recs))
                        if kept is not None:
                            pooled_strict.append((cyc, kept))
            t_score += time.perf_counter() - t0

    pooled.sort(key=lambda item: item[0])
    pooled_strict.sort(key=lambda item: item[0])
    logger.info(
        f"[perms] {gene}: {len(positions)} columns, {n_hits} discovery rows in {n_cycles} cycles; "
        f"setup {t_setup:.2f}s, kernel {t_kernel:.2f}s, entries {t_entries:.2f}s, ctx load {t_ctx:.2f}s, scoring {t_score:.2f}s")
    return pooled, pooled_strict
