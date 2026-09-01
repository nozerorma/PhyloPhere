"""Adapter: a Tier 1 pipeline run -> site-recovery metrics.

Tier 1 asks one thing: given a dataset where a published study demonstrated
convergence and *externally validated* specific sites, does the pipeline recover
them and where do they rank. This module reads a completed run's output tables
and answers that against a ``*.sites.tsv`` truth set.

Run layout (``run_pepc.sh`` / ``build_project.py``, ``POSTPROC_MODE=filter``);
``<run>`` is the ``<trait>_complete`` directory::

    <run>/caastools/discovery.tab                     candidate CAAS (per pos x scheme)
    <run>/scoring/position_scores.tsv                 per Gene x Position CAAS_score + gates
    <run>/ct_disambiguation/caas_convergence_master.csv   per (pos, scheme) ASR verdict
    <run>/selection/fade/{top,bottom}/fade_sites_*.csv    FADE BF per site (often empty)

POSITION CONVENTION.  The pipeline reports **0-based alignment column indices**
(``alimport.py``: ``range(0, alignment_length)``; ``scoring_compute.R`` even does
``site - 1`` on FADE's 1-based output to line it up).  Truth sets are **1-based**
reference coordinates.  When the fixture alignment column == reference position
(PEPC: maize P04711, 1:1) the map is ``pipeline_pos = ref_pos - 1``.  With a
reference row carried in the alignment, map ref->column with
``truthset.map_site_to_alignment`` first, then subtract 1 (pass ``ref_row``).

Residue verification.  A CAAS string is ``FG/BG`` (e.g. ``SSSS/AAAA``).  A truth
site counts as recovered only if a CAAS at the mapped position has fg mode ==
``alt_aa`` and bg mode == ``ref_aa`` (when both are known).  This guards against
truth-table transcription being off by a residue: with ``slop > 0`` the adapter
searches +/- slop columns but still requires the residues to match.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .truthset import Site, load_sites, map_site_to_alignment


# --------------------------------------------------------------------------- #
# run reader
# --------------------------------------------------------------------------- #
def _read_tsv(path: Path, delim: str = "\t") -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def _num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def _mode(seq: str) -> str:
    seq = seq.replace("-", "").replace(".", "")
    return Counter(seq).most_common(1)[0][0] if seq else ""


@dataclass
class DiscHit:
    pos0: int                       # 0-based alignment column
    schemes: set[str] = field(default_factory=set)
    caas_strings: set[str] = field(default_factory=set)
    patterns: set[str] = field(default_factory=set)
    pvalue_min: float = math.nan

    @property
    def fg_mode(self) -> str:
        return _mode("".join(s.split("/")[0] for s in self.caas_strings))

    @property
    def bg_mode(self) -> str:
        return _mode("".join(s.split("/")[1] for s in self.caas_strings if "/" in s))


class Tier1Run:
    """Lazy reader over one ``<trait>_complete`` output directory."""

    def __init__(self, run_dir: str | Path):
        self.dir = Path(run_dir)
        if not (self.dir / "caastools" / "discovery.tab").exists():
            # tolerate being handed the parent (.../out) — pick the _complete dir
            hits = sorted(self.dir.glob("**/caastools/discovery.tab"))
            if hits:
                self.dir = hits[0].parent.parent
        self._disc: dict[tuple[str, int], DiscHit] | None = None
        self._pos: dict[tuple[str, int], dict[str, str]] | None = None
        self._master: dict[tuple[str, int], list[dict[str, str]]] | None = None

    # -- discovery.tab ------------------------------------------------------- #
    def discovery(self) -> dict[tuple[str, int], DiscHit]:
        if self._disc is None:
            self._disc = {}
            for r in _read_tsv(self.dir / "caastools" / "discovery.tab"):
                key = (r["gene"], int(r["position"]))
                h = self._disc.setdefault(key, DiscHit(pos0=int(r["position"])))
                h.schemes.add(r.get("caap_group", r.get("mode", "?")))
                h.caas_strings.add(r["caas"])
                h.patterns.add(r.get("pattern", ""))
                p = _num(r.get("pvalue"))
                h.pvalue_min = p if math.isnan(h.pvalue_min) else min(h.pvalue_min, p)
        return self._disc

    # -- scoring/position_scores.tsv -------------------------------------- #
    def position_scores(self) -> dict[tuple[str, int], dict[str, str]]:
        if self._pos is None:
            p = self.dir / "scoring" / "position_scores.tsv"
            self._pos = {(r["Gene"], int(r["Position"])): r
                         for r in _read_tsv(p)} if p.exists() else {}
        return self._pos

    def scored_positions(self, gene: str | None = None) -> list[float]:
        return [_num(r["CAAS_score"]) for (g, _), r in self.position_scores().items()
                if gene is None or g == gene]

    def rank_percentile(self, gene: str, pos0: int) -> float:
        """1.0 = top-scoring position, 0.0 = lowest. Ties get the *average* rank
        (so a cluster of k positions tied at the max sit at 1 - (k-1)/2/(n-1),
        not all at 1.0)."""
        rows = self.position_scores()
        me = rows.get((gene, pos0))
        if me is None:
            return math.nan
        scores = sorted((_num(r["CAAS_score"]) for r in rows.values()), reverse=True)
        n = len(scores)
        if n <= 1:
            return 1.0
        v = _num(me["CAAS_score"])
        lo = scores.index(v)                       # first (best) 0-based rank
        hi = len(scores) - 1 - scores[::-1].index(v)  # last 0-based rank
        avg_rank = (lo + hi) / 2
        return 1.0 - avg_rank / (n - 1)

    # -- ct_disambiguation/caas_convergence_master.csv --------------------- #
    def master(self) -> dict[tuple[str, int], list[dict[str, str]]]:
        if self._master is None:
            self._master = {}
            p = self.dir / "ct_disambiguation" / "caas_convergence_master.csv"
            if p.exists():
                for r in _read_tsv(p, delim=","):
                    self._master.setdefault((r["gene"], int(r["msa_pos"])), []).append(r)
        return self._master

    # -- selection/fade/{dir}/fade_sites_{dir}.csv ------------------------- #
    def fade_sites(self, direction: str) -> dict[tuple[str, int], float]:
        p = self.dir / "selection" / "fade" / direction / f"fade_sites_{direction}.csv"
        out: dict[tuple[str, int], float] = {}
        if not p.exists():
            return out
        for r in _read_tsv(p, delim=","):
            g = r.get("gene", r.get("Gene", ""))
            raw = r.get("position", r.get("site", ""))
            if raw == "":
                continue
            # 'site' is HyPhy 1-based; 'position' is already the 0-based CAAS frame
            pos0 = int(raw) - 1 if "site" in r and "position" not in r else int(raw)
            bf = _num(r.get("max_bf", r.get("max_site_bf", r.get("bf"))))
            out[(g, pos0)] = bf
        return out


# --------------------------------------------------------------------------- #
# per-site scoring
# --------------------------------------------------------------------------- #
@dataclass
class SiteResult:
    gene: str
    ref_seq: str
    ref_pos: int                 # 1-based, truth-set coordinate
    ref_aa: str
    alt_aa: str
    tier: str
    source: str
    pipeline_pos: int            # 0-based column the truth site maps to
    recovered: bool
    residue_match: bool          # fg/bg residues matched (not just position)
    matched_pos: int | None      # actual 0-based pos of the matching CAAS (== pipeline_pos + offset)
    caas_score: float
    rank_percentile: float
    gate_sig: bool
    n_schemes: int
    scheme_set: str
    convergence_types: list[str]
    asr_path_score: float
    fade_bf: float
    in_comparator: dict[str, bool]

    def as_dict(self) -> dict:
        d = dict(self.__dict__)
        for k in ("caas_score", "rank_percentile", "asr_path_score", "fade_bf"):
            if isinstance(d[k], float) and math.isnan(d[k]):
                d[k] = None
        return d


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion — sane at small n / k∈{0,n}."""
    if n == 0:
        return (math.nan, math.nan)
    phat = k / n
    denom = 1 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def _resolve_pipeline_pos(site: Site, ref_row: str | None) -> int:
    """ref-coordinate 1-based -> pipeline 0-based column.

    ref_row given  : map through the gapped reference row, then -1 (0-based).
    ref_row None   : fixture alignment column == ref position (1:1), so -1.
    """
    if ref_row:
        return map_site_to_alignment(site.position, ref_row) - 1
    return site.position - 1


def _score_one(site: Site, run: Tier1Run, ref_row: str | None, slop: int,
               comparators: dict[str, set[int]]) -> SiteResult:
    base = _resolve_pipeline_pos(site, ref_row)
    disc = run.discovery()
    ps = run.position_scores()
    know_res = site.ref_aa not in (".", "") and site.alt_aa not in (".", "")

    matched = None
    residue_match = False
    # search base, then outward by +/- 1 .. slop.
    #   residues known -> require fg mode == alt_aa AND bg mode == ref_aa
    #   residues "."   -> positional hit at base only (slop needs residues to anchor it)
    order = [base] + [base + d for k in range(1, slop + 1) for d in (-k, k)]
    for cand in order:
        h = disc.get((site.gene, cand))
        if h is None:
            continue
        if know_res:
            if h.fg_mode == site.alt_aa and h.bg_mode == site.ref_aa:
                matched, residue_match = cand, True
                break
        elif cand == base:
            matched = cand
            break
    recovered = matched is not None

    use_pos = matched if matched is not None else base
    row = ps.get((site.gene, use_pos), {}) if matched is not None else {}
    mrows = run.master().get((site.gene, use_pos), []) if matched is not None else []
    fade = {**run.fade_sites("top"), **run.fade_sites("bottom")}

    return SiteResult(
        gene=site.gene, ref_seq=site.ref_seq, ref_pos=site.position,
        ref_aa=site.ref_aa, alt_aa=site.alt_aa, tier=site.tier, source=site.source,
        pipeline_pos=base,
        recovered=recovered,
        residue_match=residue_match,
        matched_pos=matched,
        caas_score=_num(row.get("CAAS_score")),
        rank_percentile=run.rank_percentile(site.gene, use_pos) if row else math.nan,
        gate_sig=str(row.get("gate_sig", "")).strip().upper() == "TRUE",
        n_schemes=int(_num(row.get("n_schemes"))) if row.get("n_schemes") else 0,
        scheme_set=row.get("scheme_set", ""),
        convergence_types=sorted({m.get("convergence_type", "") for m in mrows} - {""}),
        asr_path_score=_num(row.get("asr_score")),
        fade_bf=(fade.get((site.gene, use_pos), math.nan)
                 if matched is not None else math.nan),
        in_comparator={m: site.position in s for m, s in comparators.items()},
    )


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
@dataclass
class Tier1Report:
    run_dir: str
    truth_path: str
    n_truth: int
    n_scored_positions: int
    per_site: list[SiteResult]
    recall_by_tier: dict[str, dict]
    median_rank_percentile: dict[str, float]     # over recovered sites, per tier + "all"
    comparator: dict[str, dict]                  # method -> {truth_overlap, we_also_recovered, ...}
    extra_calls: list[dict]                      # discovered CAAS positions not in the truth set

    def to_dict(self) -> dict:
        def _clean(o):
            if isinstance(o, float) and math.isnan(o):
                return None
            if isinstance(o, dict):
                return {k: _clean(v) for k, v in o.items()}
            if isinstance(o, list):
                return [_clean(v) for v in o]
            return o
        return _clean({
            "run_dir": self.run_dir,
            "truth_path": self.truth_path,
            "n_truth": self.n_truth,
            "n_scored_positions": self.n_scored_positions,
            "recall_by_tier": self.recall_by_tier,
            "median_rank_percentile": self.median_rank_percentile,
            "comparator": self.comparator,
            "extra_calls": self.extra_calls,
            "per_site": [s.as_dict() for s in self.per_site],
        })


# ConDor's 5 true-positive PEPC C4 sites (Morel et al. 2024, GBE, evae040):
# M749T P540T E572Q A780S H665N — bovine... no, maize P04711 numbering.
# PCOC calls for the same fixture: TODO, transcribe from Morel 2024 SI.
DEFAULT_COMPARATORS: dict[str, dict[str, set[int]]] = {
    "pepc_c4": {"ConDor": {540, 572, 665, 749, 780}},
}


def score_tier1(
    truth_path: str | Path,
    run_dir: str | Path,
    *,
    ref_row: str | None = None,
    slop: int = 0,
    comparators: dict[str, set[int]] | None = None,
    dataset: str | None = None,
) -> Tier1Report:
    sites = load_sites(truth_path)
    run = Tier1Run(run_dir)
    comps = comparators if comparators is not None else \
        DEFAULT_COMPARATORS.get(dataset or "", {})

    per_site = [_score_one(s, run, ref_row, slop, comps) for s in sites]

    tiers = sorted({s.tier for s in per_site})
    recall_by_tier: dict[str, dict] = {}
    median_rp: dict[str, float] = {}
    for tier in tiers + ["all"]:
        grp = per_site if tier == "all" else [s for s in per_site if s.tier == tier]
        n = len(grp)
        k = sum(s.recovered for s in grp)
        lo, hi = _wilson(k, n)
        recall_by_tier[tier] = {
            "n": n, "recovered": k, "recall": (k / n if n else math.nan),
            "wilson_95ci": [round(lo, 3), round(hi, 3)],
        }
        rps = sorted(s.rank_percentile for s in grp
                     if s.recovered and not math.isnan(s.rank_percentile))
        median_rp[tier] = rps[len(rps) // 2] if rps else math.nan

    # comparator overlap
    recovered_pos = {s.ref_pos for s in per_site if s.recovered}
    truth_pos = {s.position for s in sites}
    comp_report: dict[str, dict] = {}
    for method, called in comps.items():
        comp_report[method] = {
            "n_calls": len(called),
            "in_truth": sorted(called & truth_pos),
            "we_also_recovered": sorted(called & recovered_pos),
            "they_called_we_missed": sorted((called & truth_pos) - recovered_pos),
        }

    # extra CAAS the pipeline called that are not truth sites (mapped back to ref coords)
    truth_pipe = {(s.gene, _resolve_pipeline_pos(s, ref_row)) for s in sites}
    ps = run.position_scores()
    extra = []
    for (gene, pos0), h in sorted(run.discovery().items()):
        if (gene, pos0) in truth_pipe:
            continue
        row = ps.get((gene, pos0), {})
        extra.append({
            "gene": gene, "pipeline_pos": pos0, "ref_pos_est": pos0 + 1,
            "caas": sorted(h.caas_strings)[0], "schemes": sorted(h.schemes),
            "caas_score": _num(row.get("CAAS_score")) if row else None,
            "gate_sig": str(row.get("gate_sig", "")).strip().upper() == "TRUE",
        })

    return Tier1Report(
        run_dir=str(run.dir), truth_path=str(truth_path),
        n_truth=len(sites), n_scored_positions=len(run.position_scores()),
        per_site=per_site, recall_by_tier=recall_by_tier,
        median_rank_percentile=median_rp, comparator=comp_report, extra_calls=extra,
    )
