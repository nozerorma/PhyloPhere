#!/usr/bin/env python3
"""
The permulation null per gene: replay of the labelings over each gene's ASR, with the passes that follow it.

Pass A (`process_all_genes_perms`) loads each gene's ASR once, scores every labeling and writes one detail shard per
gene; pass B (`_finalize_perm_scores`) scores each cycle against genome-wide pools and writes the null tables. Also
the readers of the detail shards and the per-gene record converter the master rows are built from.
"""

import functools
import gzip
import logging
import multiprocessing as mp
import os
import random
import sys
import time
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set, Any

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "src"))

from src.utils.concurrency import plan_concurrency, init_worker
from src.utils.amino import normalize_amino_list
from src.core.driver import load_gene_context, pool_labelings, score_labelings
from src.core.scores import DIRECTIONS, collapse_sides, direction_values, gene_scores, position_score
from src.data.models import CAASPosition

logger = logging.getLogger(__name__)


def convert_convergence_result_to_dict(
    result,
    multi_hypothesis: Optional[str],
    alignment=None,
    seq_by_id: Optional[Dict] = None,
    seq_by_species: Optional[Dict] = None,
    trait_pairs: Optional[Dict[int, List[Tuple[str, str]]]] = None,
    taxid_to_species: Optional[Dict] = None,
) -> Dict:
    """
    Convert a ConvergenceResult-like object to a JSON-serializable dict.
    This version is strictly attribute-safe: it never assumes dict APIs
    once conversion begins.

    Stability is assessed from the metadata-provided amino-encoded pattern
    (`amino_encoded`), not from reconstructed multi-caas multisets.
    """
    from types import SimpleNamespace

    # Normalize dict -> namespace for consistent attribute access
    if isinstance(result, dict):
        try:
            ns = SimpleNamespace(**result)

            # Common key mappings
            if not hasattr(ns, "position") and "position" in result:
                ns.position = result.get("position")

            if not hasattr(ns, "is_significant") and "is_significant" in result:
                ns.is_significant = result.get("is_significant")

            if not hasattr(ns, "gene") and "gene" in result:
                ns.gene = result.get("gene")

            if not hasattr(ns, "tag") and "tag" in result:
                ns.tag = result.get("tag")

            if not hasattr(ns, "caas") and "caas" in result:
                ns.caas = result.get("caas")

            if not hasattr(ns, "pair_details") and "pairs" in result:
                ns.pair_details = result.get("pairs")

            result = ns
        except Exception:
            # If normalization fails, keep original; downstream getattr will be defensive
            pass

    # Core identity
    result_dict: Dict[str, Any] = {
        "gene": getattr(result, "gene", None),
        "msa_pos": getattr(result, "msa_pos", None) or getattr(result, "position", None),  # 0-based
        "position": getattr(result, "position", None),
        "tag": getattr(result, "tag", None),
        "caas": getattr(result, "caas", None),
        "caap_group": getattr(result, "caap_group", "US"),
        "amino_encoded": getattr(result, "amino_encoded", ""),
        "multi_hypothesis": multi_hypothesis,
        # Hypotheses that drove >= 1 changed domain on THIS SIDE -- the sole
        # hypothesis-provenance column downstream reads (SCORING's pos_scores
        # aggregation consumes it by name); never nulled by a multi-hypothesis
        # collision, unlike the retired `hypothesis`-derived `trait` column.
        "participating_hypotheses": getattr(result, "participating_hypotheses", None) or "",
        # Harvest size (M) for this (position, scheme) pool.
        "n_hypotheses": getattr(result, "n_hypotheses", None),
        # Cross-hypothesis support tallies for the arbitrary first-row
        # passthrough fields above (tag/caas/amino_encoded stay as-is).
        "tag_support": getattr(result, "tag_support", "") or "",
        "caas_support": getattr(result, "caas_support", "") or "",
        "amino_encoded_support": getattr(result, "amino_encoded_support", "") or "",
    }

    # Pattern classification
    result_dict["convergence_type"] = getattr(result, "convergence_type", None)

    # First-class direction key (top / bottom / none).
    result_dict["side"] = getattr(result, "side", "none")

    # CAAS convergence score (core v3): pooled per-side domain mean;
    # derived_agreement is the diagnostic agree_num/agree_den.
    result_dict["asr_path_score"] = getattr(result, "asr_path_score", None)
    result_dict["derived_agreement"] = getattr(result, "derived_agreement", None)
    result_dict["agreement_ambiguous"] = getattr(result, "agreement_ambiguous", None)

    # ── Per-domain flat block (scoring_v2 core v3) ────────────────────────────
    # domain_<d>_score from domain_scores; domain_<d>_anc_aa / _top_aa / _bot_aa
    # from the modal harvest residues. The FOP harvest-wide, per-scheme
    # derived_agreement rebuild (V3-3/V3-4 null side) reads these.
    domain_meta = getattr(result, "domain_meta", None) or {}
    domain_scores = getattr(result, "domain_scores", None) or {}
    anc_aa = getattr(result, "domain_anc_aa", None) or {}
    der_top_aa = getattr(result, "domain_der_top_aa", None) or {}
    der_bot_aa = getattr(result, "domain_der_bot_aa", None) or {}
    der_support_top_aa = getattr(result, "domain_der_support_top_aa", None) or {}
    der_support_bot_aa = getattr(result, "domain_der_support_bot_aa", None) or {}
    anc_support_aa = getattr(result, "domain_anc_support_aa", None) or {}
    if domain_meta or domain_scores or anc_aa or der_top_aa or der_bot_aa:
        for d, meta in (domain_meta.items() if isinstance(domain_meta, dict) else []):
            m = meta or {}
            result_dict[f"domain_{d}_posterior"] = m.get("posterior")
        for d, s in domain_scores.items():
            result_dict[f"domain_{d}_score"] = s
        for d in set(anc_aa) | set(der_top_aa) | set(der_bot_aa):
            result_dict[f"domain_{d}_anc_aa"] = anc_aa.get(d)
            result_dict[f"domain_{d}_top_aa"] = der_top_aa.get(d, "")
            result_dict[f"domain_{d}_bot_aa"] = der_bot_aa.get(d, "")
        for d in set(anc_support_aa) | set(der_support_top_aa) | set(der_support_bot_aa):
            result_dict[f"domain_{d}_anc_aa_support"] = anc_support_aa.get(d, "")
            result_dict[f"domain_{d}_top_aa_support"] = der_support_top_aa.get(d, "")
            result_dict[f"domain_{d}_bot_aa_support"] = der_support_bot_aa.get(d, "")
    else:
        # Round-trip (reloaded from the aggregation DB): carry the flat
        # domain_<d>_* keys already on the input.
        src = vars(result) if hasattr(result, "__dict__") else {}
        for k, v in src.items():
            if isinstance(k, str) and k.startswith("domain_") and (
                k.endswith("_posterior")
                or k.endswith("_score") or k.endswith("_anc_aa")
                or k.endswith("_top_aa") or k.endswith("_bot_aa")
                or k.endswith("_anc_aa_support") or k.endswith("_top_aa_support")
                or k.endswith("_bot_aa_support")
            ):
                result_dict[k] = v

    return result_dict


# ═════════════════════════════════════════════════════════════════════════════
# Null-mode (permulation-excess) processing
# ═════════════════════════════════════════════════════════════════════════════
# Give CAAS a genome-wide *excess* null for FCS pathway enrichment: replay N
# permuted phenotype labelings through the FULL position pool and score each on
# the SAME ASR posteriors as the real run. ASR posteriors are phenotype-invariant
# (load_precomputed_asr is a pure function of alignment+tree), so we load them ONCE
# per gene and replay N labelings over the cached object. The scoring itself goes
# through analyze_gene_disambiguation (hence compute_asr_path_score) VERBATIM —
# observed and null share the code path, so the null is calibrated by construction.
# See docs/CAAS_PERMULATION_EXCESS.md.


def _read_resample_labelings(resample_dir: str) -> Dict[str, Tuple[List[str], List[str]]]:
    """Every resample labeling -> (fg_species, bg_species) by cycle tag (core.labelings.read_labelings)."""
    from src.core.labelings import read_labelings

    return {t: (list(l.fg), list(l.bg)) for t, l in read_labelings(resample_dir).items()}


def _read_fop_pairs(path: str) -> Dict[str, Dict[Tuple[str, int], float]]:
    """fop_pairs.tsv -> {base_cycle: {(H<n>, domain): pss_score}} (core.labelings.read_pss)."""
    from src.core.labelings import read_pss

    return read_pss(path)


def _parse_discovery_entries(
    handle,
    cycle_tags,
    gene_filter: Optional[str] = None,
) -> Dict[str, List[CAASPosition]]:
    """Parse a perm-replay-discovery TSV stream into ``{cycle_tag: [CAASPosition, ...]}``.

    Shared by ``_perms_worker_replay``'s two input layouts:

      * the single concatenated ``perm_discovery_file`` -- carries a ``gene``
        column; rows are filtered to ``gene_filter``;
      * a per-gene shard ``perm_discovery/<gene>.tsv`` -- no ``gene`` column,
        every row belongs to this gene (``gene_filter=None``).

    Both layouts otherwise share the exact same columns and the same
    ``CAASPosition`` construction, so keeping one parser here stops the two
    call-site copies from drifting as the discovery schema changes.
    """
    out: Dict[str, List[CAASPosition]] = {}
    header = handle.readline()
    if not header:
        return out
    cols = header.rstrip("\n").split("\t")
    col_indices = {col: idx for idx, col in enumerate(cols)}
    if "cycle" not in col_indices or "position" not in col_indices:
        return out
    gene_idx = col_indices.get("gene")
    if gene_filter is not None and gene_idx is None:
        return out

    cyc_idx = col_indices["cycle"]
    pos_idx = col_indices["position"]
    caas_idx = col_indices.get("caas")
    ae_idx = col_indices.get("amino_encoded")
    icm_idx = col_indices.get("is_conserved_meta")
    cp_idx = col_indices.get("conserved_pair")
    grp_idx = col_indices.get("caap_group")
    guard = max(
        i for i in (cyc_idx, pos_idx, gene_idx if gene_filter is not None else None)
        if i is not None
    )

    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= guard:
            continue
        if gene_filter is not None and parts[gene_idx].strip() != gene_filter:
            continue
        cyc = parts[cyc_idx].strip()
        if not cyc or cyc not in cycle_tags:
            continue
        try:
            pos0 = int(parts[pos_idx].strip())
        except ValueError:
            continue

        caas = parts[caas_idx] if caas_idx is not None and caas_idx < len(parts) else ""
        parts_caas = caas.split("/") if caas else []
        trait1 = normalize_amino_list(list(parts_caas[0])) if len(parts_caas) == 2 else []
        trait0 = normalize_amino_list(list(parts_caas[1])) if len(parts_caas) == 2 else []

        cp = parts[cp_idx].strip() if cp_idx is not None and cp_idx < len(parts) else ""
        if cp and ":" in cp.split(",")[0]:
            cp = cp.split(":", 1)[1]

        entry = CAASPosition(
            position=pos0,
            position_one_based=pos0 + 1,
            tag=f"POS{pos0}",
            caas=caas,
            trait1_aa=trait1,
            trait0_aa=trait0,
            caap_group=parts[grp_idx] if grp_idx is not None and grp_idx < len(parts) else "US",
            amino_encoded=parts[ae_idx] if ae_idx is not None and ae_idx < len(parts) else "",
            is_conserved_meta=parts[icm_idx] in ("TRUE", "True", "1") if icm_idx is not None and icm_idx < len(parts) else False,
            conserved_pair=cp,
        )
        out.setdefault(cyc, []).append(entry)
    return out


def build_cycle_inputs(
    perm_discovery_path: str,
    resample_dir: str,
    cycles: Optional[List[str]] = None,
) -> Tuple[List[str], Dict[str, Tuple[List[str], List[str]]]]:
    """Resolve the (fg, bg) labeling for each cycle to replay.

    Returns the cycles in-memory, straight from `_read_resample_labelings`: the
    (fg, bg) lists are already available here, so no per-cycle trait file is
    written to disk.
    """
    labelings = _read_resample_labelings(resample_dir)
    target_cycles = cycles if cycles else sorted(labelings.keys())
    cycle_labelings = {c: labelings[c] for c in target_cycles if c in labelings}

    logger.info(f"[perms] prepared trait inputs for {len(target_cycles)} cycles")
    return target_cycles, cycle_labelings


# No per-scheme weight any more. scoring_compute.R section 2g aggregates a
# position's schemes with a MEAN of caas_row, not a 0.2-weighted sum, because the
# number of detecting schemes is a biochemical-distance property of the
# substitution rather than evidence strength. The null mirrors that exactly.


@functools.lru_cache(maxsize=8)
def _scan_perm_discovery_dir(disc_path: Path) -> Dict[str, Path]:
    """One-time directory listing of a per-gene perm-discovery shard directory,
    memoized per directory. `_perms_worker_replay`'s directory-mode branch used
    to re-run `disc_path.iterdir()` (an O(N_files) NFS directory scan, ~16,100
    genes in production) on EVERY call -- once per gene pre-Stage-2, once per
    CHUNK of a gene after it (docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md),
    multiplying the same NFS metadata-call cost `find_gene_alignment`
    (io_utils.py) had -- see `_scan_alignment_dir` there for the matching fix
    and full rationale. Caching turns N calls x O(N_files) into O(N_files)
    once + N calls x O(1) dict lookup.
    """
    by_prefix: Dict[str, Path] = {}
    for p in disc_path.iterdir():
        if p.is_file():
            prefix = p.name.split(".", 1)[0]
            by_prefix.setdefault(prefix, p)
    return by_prefix


def _chunk_gene_cycles(cycle_tags: List[str], target_chunk_size: int) -> List[List[str]]:
    """Split cycle_tags into sub-chunks of ~target_chunk_size, never splitting one
    base cycle's "<base>~H*" variants across chunks -- required because
    _perms_worker_replay's FOP domain-pooling needs all of a base cycle's
    hypothesis variants together (see docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md
    Stage 2). Order-preserving; a target_chunk_size >= the input's own base-cycle
    count returns a single chunk holding all of cycle_tags."""
    from src.convergence.fop_pool import base_cycle as _bc
    groups: Dict[str, List[str]] = {}
    for tag in cycle_tags:
        groups.setdefault(_bc(tag), []).append(tag)
    chunks: List[List[str]] = []
    cur: List[str] = []
    for tags in groups.values():
        cur.extend(tags)
        if len(cur) >= target_chunk_size:
            chunks.append(cur)
            cur = []
    if cur:
        chunks.append(cur)
    return chunks or [[]]


def _perms_worker_replay(
    gene: str,
    alignment_dir: str,
    tree_file: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    cycle_tags: List[str],
    cycle_labelings: Dict[str, Tuple[List[str], List[str]]],
    perm_discovery_file: str,
    ensembl_genes: Optional[Set[str]] = None,
    fop_pairs: Optional[Dict[str, Dict[Tuple[str, int], float]]] = None,
) -> Tuple[str, Optional[List[Tuple[str, List[Any]]]]]:
    """Phase A of a chunked gene replay (see _perms_worker_finalize for phase B).

    Loads the gene's cached ASR context and replays ONLY the given cycle_tags --
    either a gene's full cycle set (one task per gene, the pre-Stage-2 shape) or
    one base-cycle-respecting sub-chunk of it (process_all_genes_perms splits a
    large gene's replay across multiple workers to bound both wall time behind
    the single largest gene and peak per-worker memory, which used to hold one
    gene's ENTIRE multi-thousand-cycle result set at once -- see
    docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md Stage 2). FOP domain-pooling of
    "<base>~H*" variants down to one record per base cycle is safe to do per-chunk
    because a chunk boundary (via _chunk_gene_cycles) never splits one base
    cycle's variants. Returns the gene name and this chunk's pooled
    (cycle_or_base, records) pairs; _perms_worker_finalize does the cross-chunk
    whole-gene reduction (n_detected, clustering, detail rows) once every chunk
    for a gene has arrived. A chunk whose context could not be loaded, or whose replay raised, returns
    (gene, None); a chunk with nothing to report returns (gene, []).
    """
    try:
        _t_ctx0 = time.perf_counter()
        ctx = load_gene_context(
            gene, alignment_dir, tree_file, taxid_mapping_path,
            asr_model, asr_cache_dir, posterior_threshold, ensembl_genes,
        )
        _t_ctx1 = time.perf_counter()
        # Chunk-sizing diagnostic (docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md
        # Stage 2): load_gene_context's cost is paid once PER CHUNK now
        # (was once per gene, pre-Stage-2), so chunk_target_size needs this
        # number to pick a chunk size where that fixed cost stays a small
        # fraction of the chunk's real replay work. Cheap (one perf_counter call
        # per chunk); left in permanently rather than removed after the first
        # measurement, since chunk_target_size may need revisiting per-cluster
        # or per-dataset (ASR cache file size varies with gene/species count).
        logger.info(
            f"[perms] {gene}: ctx load {_t_ctx1 - _t_ctx0:.3f}s "
            f"(chunk of {len(cycle_tags)} cycle-tags)"
        )
        if ctx is None:
            # load_gene_context computes ASR on a cache miss, so ctx is None only
            # when the alignment could not be found or codeml/parse failed for
            # this one gene (already warned inside).
            return (gene, None)

        # Load the gene's perm-replay discovery output in memory once. Two layouts,
        # one shared parser (_parse_discovery_entries):
        #   * a single concatenated perm_discovery_file (has a `gene` column)
        #   * a per-gene shard  perm_discovery/<gene>.tsv  (no `gene` column)
        cycle_to_entries: Dict[str, List[CAASPosition]] = {}
        disc_path = Path(perm_discovery_file)
        if disc_path.is_file():
            with open(disc_path, "r") as f:
                cycle_to_entries = _parse_discovery_entries(f, cycle_tags, gene_filter=gene)
        else:
            gene_file = _scan_perm_discovery_dir(disc_path).get(gene)
            if gene_file is not None and gene_file.exists():
                with open(gene_file, "r") as f:
                    cycle_to_entries = _parse_discovery_entries(f, cycle_tags)

        all_cycle_results = score_labelings(
            ctx, gene, cycle_tags, cycle_labelings, cycle_to_entries, posterior_threshold)
        if not all_cycle_results:
            return (gene, [])
        return (gene, pool_labelings(all_cycle_results, fop_pairs))
    except Exception as e:
        logger.error(f"[perms] replay worker failed for {gene}: {e}", exc_info=True)
        return (gene, None)


def _merge_gene_chunks(chunk_results: List[Tuple[str, List[Any]]]) -> List[Tuple[str, List[Any]]]:
    """Concatenated chunk results of one gene, in canonical (cycle tag) order.

    Chunks come back from the worker pool as they finish (imap_unordered), so their arrival order
    varies between runs. Everything written downstream follows this order (the detail shard, the
    per-cycle CAAS table, the reservoir behind perm_pos_sample.tsv), so it is fixed here. The sort
    is stable: the records of one cycle keep the order their chunk produced, and a gene replayed in
    a single chunk is left as it was (cycle tags are replayed in sorted order).
    """
    return sorted(chunk_results, key=lambda item: item[0])


def _gene_train_columns(
    train_map_index: Optional[Dict[str, Optional[str]]],
    gene: str,
    suffix: str,
    genes_without_map: List[str],
) -> Optional[Dict[int, int]]:
    """The gene's position -> untrimmed column map, or None to keep trimmed coordinates.

    None when no MAP directory was given or the gene has no MAP file; in the second case the gene is
    appended to ``genes_without_map``. A gene with an ambiguous or inconsistent MAP raises.
    """
    if train_map_index is None:
        return None
    from src.core.columns import gene_columns
    cols = gene_columns(train_map_index, gene, suffix)
    if cols is None:
        genes_without_map.append(gene)
    return cols


def _require_consistent_chunks(gene: str, failed: int, total: int) -> None:
    """A gene is replayed in `total` chunks. Every chunk failing leaves the gene out of the null, as it is left out of the observed;
    some failing leaves a null that holds only part of the cycles while N still counts all of them, so the run stops."""
    if 0 < failed < total:
        raise RuntimeError(
            f"[perms] {gene}: {failed} of {total} chunks could not be replayed (ASR or replay failure); the null would hold only part "
            "of its cycles. See the errors above.")


def _perms_worker_replay_wrapper(args):
    return _perms_worker_replay(*args)


def _perms_worker_finalize(
    gene: str,
    all_cycle_results: List[Tuple[str, List[Any]]],
    n_cycles_total: int,
    postproc_filter: bool = False,
    clust_minlen: int = 3,
    clust_maxcaas: float = 0.7,
    train_columns: Optional[Dict[int, int]] = None,
) -> Tuple[str, List[Dict[str, Any]]]:
    """Phase B of a chunked gene replay: the true whole-gene reduction over a
    gene's merged, already FOP-pooled chunk results from _perms_worker_replay --
    n_detected (needs the gene's FULL detected-cycle set),
    the CT_POSTPROC cluster filter, and detail-row emission. The output does not
    depend on how many replay chunks fed into it. n_cycles_total is passed in rather than
    derived from `all_cycle_results` because it must be the gene's (or, under FOP
    pooling, the whole run's) full cycle-tag/base-cycle universe, not just the
    subset this gene happened to detect hits in.
    """
    if not all_cycle_results:
        return (gene, [])

    # ── 1. Detection count per (position, scheme) ──────────────────────────
    # n_detected counts, per (Position, caap_group), how many of the
    # permuted-labeling cycles independently re-detected that exact CAAS.
    # It fills the detail row's own `n_detected` column (a per-position
    # replication count kept in the shards), not a
    # position-level p-value in its own right -- p.emp (scoring_compute.R,
    # from perm_pos_cycle_caas.tsv.gz) is the sole position-level permulation
    # p downstream.
    n_detected = {}
    for cyc, biochem_results in all_cycle_results:
        for r in biochem_results:
            pos = getattr(r, "position", None)
            group = getattr(r, "caap_group", "US")
            if pos is not None:
                n_detected.setdefault((pos, group), set()).add(cyc)

    n_detected_count = {key: len(cycles_set) for key, cycles_set in n_detected.items()}

    # ── 2. Emit raw per-(cycle, position, scheme) detail ──────────────────
    # Scoring itself is deliberately NOT done here. size_adj_max calibrates
    # a gene's score against its cycle's GENOME-WIDE reference pool
    # (_build_cycle_score_pools), and this worker only ever sees one gene —
    # so that pool cannot be formed at this level. The parent finalizes it
    # in pass B (see _finalize_perm_scores) once every gene's rows have
    # been counted.
    #
    # Each record is already per-side by the time it reaches here (a "both"
    # position is two records, each with its own `side` and core_s), so the
    # detail shard carries `side` directly — no OR-across-schemes step.
    # ── CT_POSTPROC cluster trains ────────────────────────────────────────
    # Per (base cycle, caap_group) core.postproc.train_flags flags this gene's
    # detected positions. The `clust` flag is emitted per detail row (0/1): pass
    # B0 reads it for the dubious-gene test, and passes B1/B2 skip clust == 1
    # rows (_is_clustered) when the trains are removed. No-op unless
    # postproc_filter is on.
    clust_by: Dict[Tuple[str, str], set] = {}
    if postproc_filter:
        from src.core.postproc import train_flags
        pos_by_cycgrp: Dict[Tuple[str, str], set] = {}
        for cyc, biochem_results in all_cycle_results:
            for r in biochem_results:
                p = getattr(r, "position", None)
                if p is None:
                    continue
                pos_by_cycgrp.setdefault(
                    (cyc, getattr(r, "caap_group", "US")), set()).add(int(p))
        clust_by = train_flags(pos_by_cycgrp, clust_maxcaas, clust_minlen, train_columns)

    detail_rows = []
    for cyc, biochem_results in all_cycle_results:
        for r in biochem_results:
            pos = getattr(r, "position", None)
            if pos is None:
                continue
            group = getattr(r, "caap_group", "US")
            asr_val = getattr(r, "asr_path_score", 0.0)
            if asr_val is None:
                asr_val = 0.0
            side = getattr(r, "side", "none") or "none"
            # `r` is already FOP-domain-pooled (or a single-contrast record)
            # and per-side by the time we get here (a "both" position is two
            # records, each with its own core_s); `side` comes off the record
            # and is the sole direction key downstream (T4b retired ct/cb).
            row = {
                "Gene": gene,
                "cycle": cyc,
                "Position": pos,
                "caap_group": group,
                "asr_path_score": asr_val,
                "n_detected": n_detected_count.get((pos, group), 1),
                "clust": 1 if int(pos) in clust_by.get((cyc, group), ()) else 0,
                "side": side,
            }
            detail_rows.append(row)

    return (gene, detail_rows)


def _null_row_caas(row: Dict[str, Any]) -> float:
    """caas_row for one null detail row (mirror of scoring_compute.R §2f).

    T1 decision E: ``caas_row = asr_score`` (no permulation percent-rank factor, on both the observed and
    null sides). Single definition shared by BOTH finalize sub-passes.
    """
    return float(row["asr_path_score"])


def _sanitize_gene_shard(gene: str) -> str:
    """Filesystem-safe stem for a per-gene detail shard. Gene ids are already
    clean in practice (Ensembl ids / HGNC symbols); this only neutralises path
    separators so a stray one cannot escape the shard directory."""
    return gene.replace(os.sep, "__").replace("/", "__").replace("\\", "__").strip() or "_"


def _is_clustered(row: Dict[str, Any]) -> bool:
    """True for a detail row whose position lies in a cluster train (`clust` = 1)."""
    return int(row.get("clust", 0) or 0) == 1


def iter_detail_rows(detail_path: Path):
    """Yield perm_pos_detail rows (dicts) from either layout:

      * a per-gene shard directory  ``perm_pos_detail/<Gene>.tsv.gz``  (current)
      * a single concatenated       ``perm_pos_detail.tsv.gz``         (legacy /
        externally supplied inputs)

    Shards are read in sorted-filename order and one shard holds exactly one
    gene, so gene-contiguity — which ``_build_cycle_score_pools`` and
    ``_finalize_perm_scores`` pass B2 both rely on — is preserved by
    construction.
    """
    import csv as _csv

    p = Path(detail_path)
    if p.is_dir():
        shards = sorted(p.glob("*.tsv.gz"))
        if not shards:
            raise FileNotFoundError(f"[perms] no *.tsv.gz shards under {p}")
        for shard in shards:
            with gzip.open(shard, "rt", newline="") as f_in:
                yield from _csv.DictReader(f_in, delimiter="\t")
    else:
        with gzip.open(p, "rt", newline="") as f_in:
            yield from _csv.DictReader(f_in, delimiter="\t")


def _cycle_gene_removal_from_detail(
    detail_path: Path,
    gene_lengths: Dict[str, float],
    mode: str,
    iqr_multiplier: float,
    extreme_percentile: float,
) -> Set[Tuple[str, str, str]]:
    """Sub-pass B0: dubious/extreme gene removal per labeling (cycle).

    One streaming read of perm_pos_detail.tsv.gz. Per (cycle, caap_group, Gene)
    accumulate the distinct detected Position count and whether any row lies in
    a train (`clust`), then apply core.postproc.gene_removal, which calibrates
    within each (cycle, caap_group) pool.
    Returns the (cycle, caap_group, Gene) units to drop from the null pool.
    """
    from src.core.postproc import GeneUnit, gene_removal

    if (mode or "none").lower() == "none":
        return set()

    seen_pos: Dict[Tuple[str, str, str], set] = {}
    has_clust: Dict[Tuple[str, str, str], bool] = {}
    for row in iter_detail_rows(detail_path):
        key = (row["cycle"], row["caap_group"], row["Gene"])
        seen_pos.setdefault(key, set()).add(int(row["Position"]))
        if int(row.get("clust", 0) or 0):
            has_clust[key] = True

    units = (
        GeneUnit(cyc, grp, gene, len(pset), has_clust.get((cyc, grp, gene), False))
        for (cyc, grp, gene), pset in seen_pos.items()
    )
    return set(gene_removal(units, gene_lengths, mode=mode,
                            iqr_multiplier=iqr_multiplier,
                            extreme_percentile=extreme_percentile))


def write_removed_units(path: Path, removed: Set[Tuple[str, str, str]]) -> None:
    """Persist the (cycle, caap_group, Gene) units dropped by gene removal, so a
    labeling's removal can be audited."""
    import csv as _csv

    with open(path, "w", newline="") as f_rm:
        w = _csv.writer(f_rm, delimiter="\t")
        w.writerow(["cycle", "caap_group", "Gene"])
        w.writerows(sorted(removed))


def _build_cycle_score_pools(
    detail_path: Path,
    removed: Optional[Set[Tuple[str, str, str]]] = None,
    remove_clusters: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """Sub-pass B1: per-cycle, per-direction pool of position-level null scores.

    The gene-level null statistic is size_adj_max = F(max)^n, mirroring
    scoring_compute.R section 4a. F must be the ECDF of the pool the gene's
    positions were actually drawn from, so every cycle needs its OWN pool --
    exactly as the observed side calibrates against the pool it discovered.
    Without this the
    null would be calibrated against a different reference than the observed
    score and the FCS p.perm comparison would be invalid.

    Directions are kept separate because scoring_compute.R uses direction-matched
    reference pools for the _top/_bottom scores.

    Pools are stored EXACTLY (sorted array.array('d')), not as a binned
    histogram. Binning was tried at 1e-4 resolution and rejected: F is raised to
    the power n, so a small ECDF error is amplified n-fold, and because the null
    score distribution is heavily tied a single bin absorbs many distinct values.
    Measured against exact pools on real detail data that gave errors up to 145%
    relative. Exactness is affordable because a pool is PER CYCLE -- ~15k
    positions, not the ~15M in the whole file -- so all 1000 cycles together cost
    on the order of 240 MB.

    Returns {cycle: {"all"|"top"|"bottom": sorted array.array('d') of scores}}.
    """
    import array
    import csv as _csv

    _rm = removed or set()
    acc: Dict[str, Dict[str, Any]] = {}

    def _per_cycle(cyc: str) -> Dict[str, Any]:
        pc = acc.get(cyc)
        if pc is None:
            pc = {k: array.array("d") for k in ("all", "top", "bottom")}
            acc[cyc] = pc
        return pc

    def _drain_side(agg: Dict[Any, Dict[str, float]]) -> None:
        # agg keyed (cyc, pos, side) -> {caap_group: caas_row}. core.scores turns each
        # into the position's score per side, then into one entry per direction: the
        # directional pools take that side's score, the global pool ONE entry per
        # position = its best side (a "both" position is not double-counted).
        by_pos: Dict[Tuple[str, int], Dict[str, float]] = {}
        for (cyc, pos, side), schemes in agg.items():
            score = position_score(schemes)
            if score is not None:
                by_pos.setdefault((cyc, pos), {})[side] = score
        for (cyc, pos), sides in by_pos.items():
            pc = _per_cycle(cyc)
            for d, v in collapse_sides(sides).items():
                pc[d].append(v)

    current_gene: Optional[str] = None
    pos_agg: Dict[Any, Dict[str, float]] = {}
    for row in iter_detail_rows(detail_path):
        gene = row["Gene"]
        if gene != current_gene:
            if current_gene is not None:
                _drain_side(pos_agg)
            current_gene = gene
            pos_agg = {}
        if _rm and (row["cycle"], row["caap_group"], gene) in _rm:
            continue
        if remove_clusters and _is_clustered(row):
            continue
        key = (row["cycle"], int(row["Position"]), row.get("side") or "none")
        pos_agg.setdefault(key, {})[row["caap_group"]] = _null_row_caas(row)
    if current_gene is not None:
        _drain_side(pos_agg)

    return {cyc: {k: array.array("d", sorted(v)) for k, v in per_cycle.items()}
            for cyc, per_cycle in acc.items()}


def _finalize_perm_scores(
    detail_path: Path,
    output_dir: Path,
    cycle_tags: List[str],
    sample_per_cycle_group: Optional[int] = None,
    removed: Optional[Set[Tuple[str, str, str]]] = None,
    seed: int = 1998,
    remove_clusters: bool = True,
) -> None:
    """Pass B: score and aggregate to gene x cycle stats.

    Mirrors scoring_compute.R's observed pipeline term for term:

        null_row_caas   = asr                                  # T1: no phen factor
        position score  = mean(null_row_caas) over that position's schemes
        gene x cycle    = size_adj_max over the cycle's positions, per direction
                          (CAAS axis; the unwired ASR axis stays on q90)

    Runs as two sub-passes over the detail file: B1 builds each cycle's
    position-score pool (_build_cycle_score_pools), B2 scores genes against it.
    The second read is unavoidable -- size_adj_max calibrates a gene against the
    distribution of its own cycle, which is not known until every position in
    that cycle has been seen.

    The detail file is gene-contiguous (the parent writes one worker result at a
    time), so a single streaming pass can aggregate per gene without ever holding
    more than one gene's rows in memory.
    """
    import csv as _csv
    import numpy as np

    scores_path = output_dir / "gene_cycle_scores.tsv"
    sample_path = output_dir / "perm_pos_sample.tsv"
    quant_path = output_dir / "perm_pos_quantiles.tsv"
    # Per-(Gene, Position, side, cycle) position score of the null cycle (core.scores), the
    # one definition every consumer reads: scoring_compute.R p.emp / SAM and the position
    # enrichment. Empty when no scheme scored the position.
    cycle_caas_path = output_dir / "perm_pos_cycle_caas.tsv.gz"
    cycle_caas_fields = ["Gene", "Position", "side", "cycle", "caas_score", "n_schemes"]

    # Reservoir size per (cycle, scheme). Bounds both the violin sample and the
    # quantile summaries at ~K * n_cycles * 5 rows regardless of run size. The
    # full detail file stays on disk, so exact statistics remain recoverable.
    if sample_per_cycle_group is None:
        sample_per_cycle_group = int(os.environ.get("CAAS_PERMS_SAMPLE_PER_CYCLE_GROUP", "200"))

    scores_fields = ["Gene", "cycle", "global_asr", "top_asr", "bottom_asr",
                     "global_caas", "top_caas", "bottom_caas"]

    rng = random.Random(seed)
    reservoirs: Dict[Tuple[str, str], List[Tuple[str, int, float, float, float]]] = {}
    seen: Dict[Tuple[str, str], int] = {}
    _rm = removed or set()

    # ── Sub-pass B1: per-cycle reference pools for the size-adjusted max ──────
    cycle_pools = _build_cycle_score_pools(detail_path, removed=_rm, remove_clusters=remove_clusters)
    logger.info("[perms] pass B1 done: size-adjust reference pools built for %d cycles",
                len(cycle_pools))

    # The detail shard always carries a `side` column: a "both" position is two
    # rows (one per side, each its own core_s); key per (cyc, pos, side) and
    # dedup the global pool by max(side) per position (T3-doc §12).

    def _q90(vals) -> float:
        return float(np.percentile(vals, 90)) if vals else 0.0

    def _na(x) -> Any:
        return "NA" if x is None else x

    def _flush(gene: str, pos_scheme: Dict[Any, Dict[str, float]], writer) -> None:
        # pos_scheme keyed (cyc, pos, side) -> {caap_group: caas_row}. Per-cycle CAAS
        # numerator/denominator first (scoring_compute.R p.emp redoes the division).
        by_pos: Dict[Tuple[str, int], Dict[str, float]] = {}
        cc_rows = []
        for (cyc, pos, side), schemes in pos_scheme.items():
            score = position_score(schemes)
            cc_rows.append({"Gene": gene, "Position": pos, "side": side, "cycle": cyc,
                            "caas_score": score, "n_schemes": len(schemes)})
            if score is not None:
                by_pos.setdefault((cyc, pos), {})[side] = score
        if writer_cc is not None and cc_rows:
            writer_cc.writerows(cc_rows)

        positions_by_cycle: Dict[str, List[Dict[str, float]]] = {}
        for (cyc, _pos), sides in by_pos.items():
            positions_by_cycle.setdefault(cyc, []).append(sides)

        no_pool = {d: () for d in DIRECTIONS}
        rows = []
        for cyc in cycle_tags:
            positions = positions_by_cycle.get(cyc, [])
            vals = direction_values(positions)
            caas = gene_scores(positions, cycle_pools.get(cyc, no_pool))
            # The ASR axis is not wired into any ranking downstream; it stays on q90
            # (0 when empty), like the observed gene_caas_asr_* columns. asr_path_score
            # is the caas_row, so it shares the position scores.
            # The CAAS axis is what FCS consumes (caas_corStat_byrank). A cycle with no
            # scored position in a direction has no score there (NA); scoring_caas_perms.R
            # fills those cells with 0 when it builds the dense genes x cycles matrix.
            rows.append({
                "Gene": gene, "cycle": cyc,
                "global_asr": _q90(vals["all"]), "top_asr": _q90(vals["top"]), "bottom_asr": _q90(vals["bottom"]),
                "global_caas": _na(caas["all"]), "top_caas": _na(caas["top"]), "bottom_caas": _na(caas["bottom"]),
            })
        writer.writerows(rows)

    n_rows = 0
    with open(scores_path, "w", newline="") as f_scores, \
         gzip.open(cycle_caas_path, "wt", newline="") as f_cc:
        reader = iter_detail_rows(detail_path)
        writer_scores = _csv.DictWriter(f_scores, fieldnames=scores_fields, delimiter="\t")
        writer_scores.writeheader()
        writer_cc = _csv.DictWriter(f_cc, fieldnames=cycle_caas_fields, delimiter="\t")
        writer_cc.writeheader()

        current_gene: Optional[str] = None
        pos_scheme: Dict[Tuple[str, int, str], Dict[str, float]] = {}

        for row in reader:
            gene = row["Gene"]
            if gene != current_gene:
                if current_gene is not None:
                    _flush(current_gene, pos_scheme, writer_scores)
                current_gene = gene
                pos_scheme = {}

            cyc = row["cycle"]
            grp = row["caap_group"]
            if _rm and (cyc, grp, gene) in _rm:
                continue
            if remove_clusters and _is_clustered(row):
                continue
            pos = int(row["Position"])
            asr = float(row["asr_path_score"])

            rc = asr  # T1 decision E: caas_row = asr_score (no phen factor)

            skey = (cyc, pos, row.get("side") or "none")
            pos_scheme.setdefault(skey, {})[grp] = rc

            # Reservoir sample stratified by (cycle, scheme): every cycle
            # contributes up to the same K rows regardless of how many detections
            # it produced, so each cycle is equally represented in the report's
            # distribution plots and per-cycle summaries can be formed. Stratifying
            # on the cycle (rather than on the gene) is what preserves per-cycle
            # resolution; sampling per gene would weight genes by how many cycles
            # happened to hit them.
            key = (cyc, grp)
            seen[key] = seen.get(key, 0) + 1
            res = reservoirs.setdefault(key, [])
            item = (gene, pos, asr, rc)
            if len(res) < sample_per_cycle_group:
                res.append(item)
            else:
                j = rng.randrange(seen[key])
                if j < sample_per_cycle_group:
                    res[j] = item
            n_rows += 1

        if current_gene is not None:
            _flush(current_gene, pos_scheme, writer_scores)

    # ── Sample + quantile summaries ────────────────────────────────────────────
    sample_fields = ["Gene", "Position", "caap_group", "cycle",
                     "asr_path_score", "null_row_caas"]
    with open(sample_path, "w", newline="") as f_sample:
        writer_sample = _csv.DictWriter(f_sample, fieldnames=sample_fields, delimiter="\t")
        writer_sample.writeheader()
        for (cyc, grp), res in reservoirs.items():
            writer_sample.writerows({
                "Gene": g, "Position": p, "caap_group": grp, "cycle": cyc,
                "asr_path_score": a, "null_row_caas": rc,
            } for (g, p, a, rc) in res)

    quant_levels = [5, 10, 25, 50, 75, 90, 95]
    quant_fields = (["cycle", "caap_group", "metric", "n_sampled", "n_total", "mean"]
                    + [f"q{q}" for q in quant_levels])
    with open(quant_path, "w", newline="") as f_quant:
        writer_quant = _csv.DictWriter(f_quant, fieldnames=quant_fields, delimiter="\t")
        writer_quant.writeheader()
        for (cyc, grp), res in reservoirs.items():
            if not res:
                continue
            for metric, idx in (("asr_path_score", 2), ("null_row_caas", 3)):
                vals = np.asarray([r[idx] for r in res], dtype=float)
                rec = {"cycle": cyc, "caap_group": grp, "metric": metric,
                       "n_sampled": len(vals), "n_total": seen.get((cyc, grp), len(vals)),
                       "mean": float(vals.mean())}
                for q, v in zip(quant_levels, np.percentile(vals, quant_levels)):
                    rec[f"q{q}"] = float(v)
                writer_quant.writerow(rec)

    logger.info(
        "[perms] pass B done: scored %d rows -> %s; sample=%s quantiles=%s cycle_caas=%s",
        n_rows, scores_path.name, sample_path.name, quant_path.name, cycle_caas_path.name,
    )


def process_all_genes_perms(
    genes: List[str],
    alignment_dir: str,
    tree_file: str,
    perm_discovery_file: str,
    resample_dir: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    workers: Optional[int],
    output_dir: Path,
    ensembl_genes_file: Optional[str] = None,
    cycles: Optional[List[str]] = None,
    max_tasks_per_child: Optional[int] = None,
    fop_pairs_file: Optional[str] = None,
    gene_lengths_file: Optional[str] = None,
    clust_minlen: int = 3,
    clust_maxcaas: float = 0.7,
    gene_filter_mode: str = "none",
    iqr_multiplier: float = 3.0,
    extreme_percentile: float = 0.99,
    postproc_filter: bool = False,
    gene_sizes: Optional[Dict[str, int]] = None,
    chunk_threshold: Optional[int] = None,
    chunk_target_size: Optional[int] = None,
    seed: int = 1998,
    remove_clusters: bool = True,
    train_map_dir: Optional[str] = None,
    train_map_suffix: str = ".map.tsv",
    detail_only: bool = False,
) -> Path:
    """Genome-wide CAAS permulation null: load ASR once per gene, replay N permuted
    labelings, and score them the same way the observed pipeline scores itself.

    Two passes, because each gene's score is calibrated against its cycle's genome-wide pool of position
    scores while workers only ever see a single gene:

      Pass A (parallel, one gene per worker) replays the labelings and streams raw
        per-(gene, cycle, position, scheme) detail to perm_pos_detail/<Gene>.tsv.gz
        (one gz shard per gene, since each worker result is one gene's complete
        row list).
      Pass B (single process, streaming) re-reads the detail shards via
        iter_detail_rows(), builds each cycle's pool of position scores and derives
        null_row_caas -> position scores -> per-(gene, cycle) size-adjusted max.

    Outputs:
      - output_dir/gene_cycle_scores.tsv     (feeds caas_perms.rds)
      - output_dir/perm_pos_cycle_caas.tsv.gz (per (Gene, Position, side, cycle)
                                              caas_score, the core.scores position
                                              score of the cycle; the R side takes
                                              the max over sides for the pooled
                                              p.emp, the sole position-level
                                              permulation p)
      - output_dir/perm_pos_detail/<Gene>.tsv.gz  (one shard per gene; re-scoring
                                              needs no ASR replay, and re-aggregation
                                              stays at one-gene peak RAM)
      - output_dir/perm_pos_detail.manifest.tsv   (Gene, n_rows per shard)
      - output_dir/perm_pos_quantiles.tsv    (per (cycle, scheme) distribution shape)
      - output_dir/perm_pos_sample.tsv       (cycle-stratified sample for violins)

    With detail_only the run stops after pass A: only perm_pos_detail/ and its manifest are written, for
    a caller that scores the union of several runs' shards in one later pass B.
    """
    import csv as _csv

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    effective_workers, _ = plan_concurrency(workers, 1, logger)

    from src.data.loaders import load_ensembl_genes
    ensembl_genes: Optional[Set[str]] = None
    if ensembl_genes_file:
        try:
            ensembl_genes = load_ensembl_genes(Path(ensembl_genes_file)) or set()
        except Exception as exc:
            logger.warning(f"[perms] failed to load ensembl genes: {exc}")

    if ensembl_genes is not None:
        genes = [g for g in genes if g in ensembl_genes]

    # FOP mirror: when resample_fop_pairs.tsv is present the resample dir carries
    # "<base>~H<m>" hypothesis labelings (resample_fop.tab) instead of / alongside
    # the plain resample_*.tab. build_cycle_inputs resolves one labeling per
    # hypothesis tag; _perms_worker_replay then domain-pools them back to one score
    # per base cycle. Parse the per-(hypothesis, domain) PSS weights once here.
    fop_pairs: Optional[Dict[str, Dict[Tuple[str, int], float]]] = None
    if fop_pairs_file and Path(fop_pairs_file).exists():
        fop_pairs = _read_fop_pairs(fop_pairs_file)
        logger.info(f"[perms] FOP mirror ON: PSS weights for {len(fop_pairs)} base cycles")

    cycle_tags, cycle_labelings = build_cycle_inputs(
        perm_discovery_file, resample_dir, cycles
    )
    if not cycle_tags:
        raise RuntimeError("[perms] no usable cycles (check export_perm_discovery + resample dir)")
    logger.info(
        f"[perms] replaying {len(cycle_tags)} cycles over {len(genes)} genes "
        f"with {effective_workers} workers"
    )

    if max_tasks_per_child is not None:
        maxtasks = int(max_tasks_per_child)
    else:
        maxtasks = int(os.environ.get("CAAS_MAX_TASKS_PER_CHILD", "50"))

    # One gz shard per gene rather than one monolithic file: Pass A already
    # iterates gene-by-gene (imap_unordered yields one gene's complete row list
    # at a time), so this is a writer-only change. Keeps every downstream
    # re-aggregation (CAAS_CORE_MERGE, this same finalizer) at one-gene peak
    # RAM instead of streaming a single multi-GB file, and lets consumers read
    # shard-by-shard. `iter_detail_rows()` is the matching reader.
    detail_dir = output_dir / "perm_pos_detail"
    detail_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "perm_pos_detail.manifest.tsv"

    # The detail shard carries `side` (a "both" position is two detail rows,
    # one per side, each with its own core_s).
    detail_fields = ["Gene", "cycle", "Position", "caap_group", "asr_path_score",
                     "n_detected", "clust", "side"]

    # Gap B: CT_POSTPROC filtering of the null candidate pool. Off by default so
    # the non-postproc null path is byte-identical; the nextflow layer flips it
    # on (params.caas_perms_postproc) to distribution-match the observed
    # filtered_discovery.tsv -> scoring_compute.R chain.
    gene_lengths: Dict[str, float] = {}
    if postproc_filter and gene_lengths_file:
        try:
            from src.core.postproc import load_gene_lengths
            gene_lengths = load_gene_lengths(gene_lengths_file)
            logger.info(f"[perms] CT_POSTPROC filter ON: {len(gene_lengths)} gene lengths, "
                        f"cluster minlen={clust_minlen} maxcaas={clust_maxcaas}, "
                        f"gene_filter_mode={gene_filter_mode}")
        except Exception as exc:
            logger.warning(f"[perms] could not load gene lengths ({exc}); "
                           "extreme-gene filter disabled")

    # Optional untrimmed-coordinate trains: with a MAP directory, each gene's cluster trains measure their span
    # in untrimmed alignment columns (core.columns). A gene without a MAP file keeps the trimmed coordinates.
    train_map_index: Optional[Dict[str, Optional[str]]] = None
    genes_without_map: List[str] = []
    if postproc_filter and train_map_dir:
        from src.core.columns import index_files
        train_map_index = index_files(train_map_dir, train_map_suffix)
        logger.info(f"[perms] cluster trains in untrimmed coordinates: {len(train_map_index)} MAP files in {train_map_dir}")

    # n_cycles_total is the SAME for every gene: all genes replay the same global
    # cycle pool (cycle_tags/cycle_labelings above), so it's computed once here
    # rather than re-derived per gene/chunk in _perms_worker_finalize.
    if fop_pairs is not None:
        from src.convergence.fop_pool import base_cycle as _bc_global
        n_cycles_total = len({_bc_global(c) for c in cycle_tags})
    else:
        n_cycles_total = len(cycle_tags)

    if chunk_threshold is None:
        chunk_threshold = int(os.environ.get("CAAS_PERMS_CHUNK_THRESHOLD", "5000"))
    if chunk_target_size is None:
        chunk_target_size = int(os.environ.get("CAAS_PERMS_CHUNK_TARGET_SIZE", "1000"))

    # Stage 2 (docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md): split a gene's
    # replay into base-cycle-respecting sub-chunks, each dispatched as its own
    # worker task, when its LPT size proxy (gene_sizes, the row-count tally
    # disambiguation_perms_main.py's _genes_from_export already computes for LPT
    # ordering) exceeds chunk_threshold. This bounds both wall time behind a
    # single oversized gene AND per-worker memory, which previously held one
    # gene's entire multi-thousand-cycle result set at once regardless of size
    # (the OOM this Stage was built to fix). Below threshold (or with no
    # gene_sizes given), a gene gets exactly one chunk -- its full cycle_tags --
    # matching the pre-Stage-2 one-task-per-gene dispatch exactly. All genes
    # replay the same global cycle_tags, so the split points are identical
    # across every chunked gene and computed once, not per gene.
    big_chunks = None
    if gene_sizes and any(gene_sizes.get(g, 0) > chunk_threshold for g in genes):
        big_chunks = _chunk_gene_cycles(cycle_tags, chunk_target_size)

    tasks: List[Tuple[str, List[str], float]] = []  # (gene, chunk_cycle_tags, est_workload)
    chunks_per_gene: Dict[str, int] = {}
    for gene in genes:
        size = (gene_sizes or {}).get(gene, 0)
        use_chunks = big_chunks if (big_chunks and size > chunk_threshold) else [cycle_tags]
        chunks_per_gene[gene] = len(use_chunks)
        for chunk in use_chunks:
            tasks.append((gene, chunk, size / len(use_chunks)))

    # LPT: dispatch the largest-estimated-workload chunk first. `genes` already
    # arrives LPT-ordered (largest gene first) from _genes_from_export; sorting
    # the flattened per-chunk task list preserves that property when a gene has
    # been split into several chunks of roughly equal size.
    tasks.sort(key=lambda t: -t[2])

    args_generator = (
        (gene, alignment_dir, tree_file, taxid_mapping_path, asr_model,
         asr_cache_dir, posterior_threshold,
         chunk, cycle_labelings, perm_discovery_file, None, fop_pairs)
        for gene, chunk, _est in tasks
    )

    # ── Pass A: replay labelings and stream the per-gene detail shards ──────────
    cycles_seen: Set[str] = set()
    pool = mp.Pool(
        processes=effective_workers, maxtasksperchild=maxtasks,
        initializer=init_worker, initargs=(1, None),
    )
    n_genes = 0
    n_detail_rows = 0
    manifest_rows: List[Tuple[str, int]] = []
    # Chunk bookkeeping: a gene's pooled per-chunk results accumulate here until
    # every one of its chunks (chunks_per_gene[gene]) has arrived, at which point
    # _perms_worker_finalize runs once and the entry is dropped -- so at most a
    # handful of genes (bounded by effective_workers, since that's how many
    # chunks can be in flight at once) are ever partially resident here.
    pending_pooled: Dict[str, List[Tuple[str, List[Any]]]] = {}
    received: Dict[str, int] = {}
    failed_chunks: Dict[str, int] = {}
    try:
        # chunksize=1: each item here is one gene-cycle-chunk's replay (real
        # seconds of work), not the many-cheap-tasks shape chunksize>1 is for.
        # A chunksize >= len(tasks) bundles ALL tasks into ONE chunk handed to a
        # SINGLE worker — the other `effective_workers - 1` workers never receive
        # any work at all for the whole batch. Confirmed via real-data timing
        # (docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md, multi-worker contention
        # investigation): with chunksize=10 and 3 genes queued, only 1 of 3
        # workers ever ran (only one gene's "TAXONOMY CONFLICT" log line fired,
        # instead of one per worker) and wall time was ~1.7x worse than
        # chunksize=1's genuinely-parallel 3-worker run. In production, 20 genes
        # / chunksize=10 = exactly 2 chunks, so only 2 of the 8 allocated workers
        # ever ran per batch regardless of pool size — mechanically capping
        # utilization at 25% before any other inefficiency, matching the ~1-2%
        # aggregate CPU utilization measured on live production jobs. This was a
        # regression from bf92df7 (2026-07-16), which replaced a per-gene
        # apply_async dispatch (no chunking, so no such cap existed) with
        # imap_unordered.
        results_iterator = pool.imap_unordered(_perms_worker_replay_wrapper, args_generator, chunksize=1)

        for _gene, chunk_pooled in results_iterator:
            if chunk_pooled is None:
                failed_chunks[_gene] = failed_chunks.get(_gene, 0) + 1
                chunk_pooled = []
            pending_pooled.setdefault(_gene, []).extend(chunk_pooled)
            received[_gene] = received.get(_gene, 0) + 1
            if received[_gene] < chunks_per_gene.get(_gene, 1):
                continue  # more chunks still in flight for this gene
            _require_consistent_chunks(_gene, failed_chunks.pop(_gene, 0), chunks_per_gene.get(_gene, 1))

            gene_pooled = _merge_gene_chunks(pending_pooled.pop(_gene))
            received.pop(_gene, None)
            _gene, detail_rows = _perms_worker_finalize(
                _gene, gene_pooled, n_cycles_total,
                postproc_filter, clust_minlen, clust_maxcaas,
                _gene_train_columns(train_map_index, _gene, train_map_suffix, genes_without_map),
            )
            if not detail_rows:
                continue

            shard_path = detail_dir / f"{_sanitize_gene_shard(_gene)}.tsv.gz"
            with gzip.open(shard_path, "wt", newline="") as f_detail:
                writer_detail = _csv.DictWriter(f_detail, fieldnames=detail_fields, delimiter="\t")
                writer_detail.writeheader()
                writer_detail.writerows(detail_rows)
            manifest_rows.append((_gene, len(detail_rows)))

            n_detail_rows += len(detail_rows)
            cycles_seen.update(row["cycle"] for row in detail_rows)
            n_genes += 1
    finally:
        pool.close()
        pool.join()

    with open(manifest_path, "w", newline="") as f_man:
        w = _csv.writer(f_man, delimiter="\t")
        w.writerow(["Gene", "n_rows"])
        w.writerows(sorted(manifest_rows))

    logger.info(
        f"[perms] pass A done: {n_genes} genes, {n_detail_rows} (gene,cycle,position,scheme) "
        f"rows across {len(cycles_seen)} cycles -> {detail_dir.name}/ ({n_genes} shards)"
    )
    if n_genes < len(genes):
        logger.info(
            f"[perms] pass A scored {n_genes}/{len(genes)} genes with per-cycle CAAS — "
            f"{len(genes) - n_genes} contributed nothing (no CAAS survived any replayed "
            f"labeling, or the alignment / ASR failed for that gene)."
        )
    if train_map_index is not None and genes_without_map:
        logger.warning(f"[perms] {len(genes_without_map)} genes have no MAP file and keep trimmed-coordinate "
                       f"trains, e.g. {sorted(genes_without_map)[:5]}")
    if n_detail_rows == 0:
        logger.error(
            "[perms] pass A produced ZERO detail rows — the permulation null is empty. "
            "Downstream caas_perms.rds / FCS p.perm will be degenerate. Check the "
            "per-cycle perm-replay discovery (export_perm_discovery) and the ASR cache."
        )

    if detail_only:
        return output_dir

    # ── Pass B: score each cycle against its own pool, aggregate ────────────────

    # ── Sub-pass B0 (Gap B): cycle-aware dubious/extreme gene removal ──────────
    removed: Set[Tuple[str, str, str]] = set()
    if postproc_filter and (gene_filter_mode or "none").lower() != "none":
        removed = _cycle_gene_removal_from_detail(
            detail_dir, gene_lengths, gene_filter_mode,
            iqr_multiplier, extreme_percentile,
        )
        logger.info("[perms] pass B0: %d (cycle, group, gene) units removed "
                    "(mode=%s) — mirrors CAAS_FILTER_GENES on the null pool",
                    len(removed), gene_filter_mode)
        write_removed_units(output_dir / "removed_units.tsv", removed)

    # _finalize_perm_scores aggregates the base-cycle-keyed detail shards. Under
    # the FOP mirror, build_cycle_inputs' `cycle_tags` are the "<base>~H<m>"
    # hypothesis-replay tags, but _perms_worker_replay domain-pools those down to ONE
    # record per base cycle before it writes any detail row (see the "FOP
    # domain-pooling" block above; it also does `cycle_tags = {_bc(c) ...}`
    # locally). Pass the base-collapsed list here too, or _flush's
    # `for cyc in cycle_tags` loop never matches the detail's `cycle` column and
    # every gene x cycle row is written as a false structural zero -> an all-zero
    # gene_cycle_scores.tsv / caas_perms.rds and a degenerate FCS p.perm.
    finalize_cycle_tags = cycle_tags
    if fop_pairs is not None:
        from src.convergence.fop_pool import base_cycle as _bc
        finalize_cycle_tags = sorted({_bc(c) for c in cycle_tags})
        logger.info("[perms] FOP mirror: collapsed %d '<base>~H*' replay tags to "
                    "%d base cycles for gene x cycle aggregation",
                    len(cycle_tags), len(finalize_cycle_tags))

    _finalize_perm_scores(
        detail_path=detail_dir,
        output_dir=output_dir,
        cycle_tags=finalize_cycle_tags,
        removed=removed,
        seed=seed,
        remove_clusters=remove_clusters,
    )

    logger.info(f"[perms] successfully aggregated {n_genes} genes to summaries inside {output_dir}")
    return output_dir


# Backward-compatible alias
