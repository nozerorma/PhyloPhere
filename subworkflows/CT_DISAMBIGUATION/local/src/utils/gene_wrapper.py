#!/usr/bin/env python3
# gene_wrapper.py — Permulation null per gene: replay of the labelings over each gene's ASR and the passes after it.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
The permulation null per gene: replay of the labelings over each gene's ASR, with the passes that follow it.

Pass A (`process_all_genes_perms`) loads each gene's ASR once, scores every labeling and writes one detail shard per
gene; pass B (`_finalize_perm_scores`) scores each cycle against genome-wide pools and writes the null tables. Also
the readers of the detail shards and the per-gene record converter the master rows are built from.

The null gives CAAS a genome-wide excess null for FCS pathway enrichment: N permuted phenotype labelings are replayed
through the full position pool and scored on the same ASR posteriors as the observed run. The posteriors depend only
on alignment and tree (not on the phenotype), so they are loaded once per gene and reused for every labeling. The
scoring goes through the same code as the observed run, so the null is calibrated against the observed by construction.

Imported by: disambiguation_perms_main.py (process_all_genes_perms), reaggregate_perm_scores.py (pass B and the
             detail readers), src/core/observed.py (convert_convergence_result_to_dict)
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
    """Convert a ConvergenceResult-like object (or a dict of its fields) to a JSON-serializable dict.

    Attribute-safe: after the dict-to-namespace normalization every field is read with getattr and a
    default, so a record that lacks a field yields None (or the stated default) instead of raising.
    """
    from types import SimpleNamespace

    # A dict input is turned into a namespace so that every field is read the same way
    if isinstance(result, dict):
        try:
            ns = SimpleNamespace(**result)

            # Keys whose attribute name differs from the dict key (pairs -> pair_details)
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
            # Keep the original object: every later read is a defensive getattr
            pass

    # Identity and pattern of the position
    result_dict: Dict[str, Any] = {
        "gene": getattr(result, "gene", None),
        "msa_pos": getattr(result, "msa_pos", None) or getattr(result, "position", None),  # 0-based
        "position": getattr(result, "position", None),
        "tag": getattr(result, "tag", None),
        "caas": getattr(result, "caas", None),
        "caap_group": getattr(result, "caap_group", "US"),
        "amino_encoded": getattr(result, "amino_encoded", ""),
        "multi_hypothesis": multi_hypothesis,
        # Hypotheses that drove at least one changed domain on THIS SIDE: the hypothesis-provenance
        # column the scoring reads (scoring_compute.R aggregates it by name). It is never blanked
        # when several hypotheses collide on the position.
        "participating_hypotheses": getattr(result, "participating_hypotheses", None) or "",
        # Number of hypotheses (M) harvested for this (position, scheme) pool.
        "n_hypotheses": getattr(result, "n_hypotheses", None),
        # Cross-hypothesis support tallies for tag, caas and amino_encoded, which above keep the
        # value of the first row.
        "tag_support": getattr(result, "tag_support", "") or "",
        "caas_support": getattr(result, "caas_support", "") or "",
        "amino_encoded_support": getattr(result, "amino_encoded_support", "") or "",
    }

    # Pattern classification
    result_dict["convergence_type"] = getattr(result, "convergence_type", None)

    # Direction of the record (top / bottom / none); the key every downstream step groups by.
    result_dict["side"] = getattr(result, "side", "none")

    # CAAS convergence score: the pooled per-side domain mean. derived_agreement is the
    # diagnostic agree_num/agree_den.
    result_dict["asr_path_score"] = getattr(result, "asr_path_score", None)
    result_dict["derived_agreement"] = getattr(result, "derived_agreement", None)
    result_dict["agreement_ambiguous"] = getattr(result, "agreement_ambiguous", None)

    # ── Per-domain flat block ─────────────────────────────────────────────────
    # domain_<d>_score from domain_scores; domain_<d>_anc_aa / _top_aa / _bot_aa (and the
    # *_support tallies) from the modal harvest residues of domain d.
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
        # The input is a record that was already flattened (e.g. read back from a table):
        # carry its domain_<d>_* keys over unchanged.
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


# ── Permulation null: readers, chunking and the per-gene replay ───────────────


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
    """Parse a perm-replay discovery TSV stream into ``{cycle_tag: [CAASPosition, ...]}``.

    Serves the two input layouts of ``_perms_worker_replay``:

      * the single concatenated ``perm_discovery_file``, which carries a ``gene``
        column; rows are filtered to ``gene_filter``;
      * a per-gene shard ``perm_discovery/<gene>.tsv``, which has no ``gene`` column:
        every row belongs to this gene (``gene_filter=None``).

    The layouts share their columns, so one parser builds the ``CAASPosition`` objects of both.
    Rows of cycles outside ``cycle_tags`` and rows with a non-integer position are skipped.
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
    """Resolve the (fg, bg) labeling of each cycle to replay.

    Returns (cycle tags, {tag: (fg, bg)}), all in memory: `cycles` when given, else every
    cycle of `resample_dir`; a requested tag missing from the resample files is dropped from
    the labelings. No per-cycle trait file is written. `perm_discovery_path` is not read here.
    """
    labelings = _read_resample_labelings(resample_dir)
    target_cycles = cycles if cycles else sorted(labelings.keys())
    cycle_labelings = {c: labelings[c] for c in target_cycles if c in labelings}

    logger.info(f"[perms] prepared trait inputs for {len(target_cycles)} cycles")
    return target_cycles, cycle_labelings


# A position's score is the mean of caas_row over its schemes (core.scores.position_score), not a
# weighted sum: the number of schemes that detect a substitution reflects its biochemical distance
# rather than the strength of the evidence. The observed side (scoring_compute.R) uses the same score.


@functools.lru_cache(maxsize=8)
def _scan_perm_discovery_dir(disc_path: Path) -> Dict[str, Path]:
    """Map the gene prefix of every file in a per-gene perm-discovery shard directory to its path.

    Memoized per directory: `_perms_worker_replay` calls it once per gene chunk, and listing the
    directory each time would repeat an O(N_files) scan on the file system. The key is the file name
    up to the first dot (the same convention as `io_utils._scan_alignment_dir`).
    """
    by_prefix: Dict[str, Path] = {}
    for p in disc_path.iterdir():
        if p.is_file():
            prefix = p.name.split(".", 1)[0]
            by_prefix.setdefault(prefix, p)
    return by_prefix


def _chunk_gene_cycles(cycle_tags: List[str], target_chunk_size: int) -> List[List[str]]:
    """Split cycle_tags into sub-chunks of about target_chunk_size tags.

    The "<base>~H*" variants of one base cycle never straddle two chunks, because the FOP
    domain-pooling of _perms_worker_replay needs all of them together. Order-preserving; a
    target_chunk_size at least as large as the number of tags returns a single chunk."""
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
    """Replay step of a chunked gene replay (_perms_worker_finalize is the reduction step).

    Loads the gene's ASR context and replays ONLY the given cycle_tags: the gene's full cycle set,
    or one base-cycle-respecting sub-chunk of it (process_all_genes_perms splits the replay of a
    large gene across workers to bound the wall time behind the largest gene and the memory a
    worker holds). The FOP domain-pooling of "<base>~H*" variants into one record per base cycle
    can be done per chunk because _chunk_gene_cycles never splits the variants of a base cycle.
    Returns the gene name and the pooled (cycle_or_base, records) pairs of this chunk;
    _perms_worker_finalize does the whole-gene reduction (n_detected, clustering, detail rows)
    once every chunk of the gene has arrived. A chunk whose context could not be loaded, or whose
    replay raised, returns (gene, None); a chunk with nothing to report returns (gene, []).
    """
    try:
        _t_ctx0 = time.perf_counter()
        ctx = load_gene_context(
            gene, alignment_dir, tree_file, taxid_mapping_path,
            asr_model, asr_cache_dir, posterior_threshold, ensembl_genes,
        )
        _t_ctx1 = time.perf_counter()
        # The context load is paid once per chunk. Logging its time shows whether chunk_target_size
        # keeps that fixed cost small against the replay work of a chunk; the ASR cache file size
        # varies with the gene and the number of species.
        logger.info(
            f"[perms] {gene}: ctx load {_t_ctx1 - _t_ctx0:.3f}s "
            f"(chunk of {len(cycle_tags)} cycle-tags)"
        )
        if ctx is None:
            # load_gene_context computes ASR on a cache miss, so ctx is None only
            # when the alignment could not be found or codeml/parse failed for
            # this one gene (already warned inside).
            return (gene, None)

        # Load the gene's perm-replay discovery rows in memory once. Two layouts, one parser
        # (_parse_discovery_entries):
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
    """Reduction step of a chunked gene replay: the whole-gene pass over the merged, already
    FOP-pooled chunk results of _perms_worker_replay.

    Computes n_detected (which needs the gene's full set of detected cycles), the CT_POSTPROC
    cluster flags and the detail rows. The output does not depend on how many replay chunks fed
    into it. The `n_cycles_total` argument is the full cycle (or, under FOP pooling, base-cycle)
    universe of the run and is not used inside this function.
    """
    if not all_cycle_results:
        return (gene, [])

    # ── 1. Detection count per (position, scheme) ─────────────────────────────
    # n_detected counts, per (Position, caap_group), how many permuted-labeling cycles
    # re-detected that CAAS. It only fills the `n_detected` column of the detail rows; the
    # position-level permulation p is p.emp of scoring_compute.R, computed from
    # perm_pos_cycle_caas.tsv.gz.
    n_detected = {}
    for cyc, biochem_results in all_cycle_results:
        for r in biochem_results:
            pos = getattr(r, "position", None)
            group = getattr(r, "caap_group", "US")
            if pos is not None:
                n_detected.setdefault((pos, group), set()).add(cyc)

    n_detected_count = {key: len(cycles_set) for key, cycles_set in n_detected.items()}

    # ── 2. Emit raw per-(cycle, position, scheme) detail ──────────────────────
    # Gene scores are not computed here. size_adj_max calibrates a gene's score against the
    # genome-wide reference pool of its cycle (_build_cycle_score_pools), and a worker sees
    # one gene only, so the pool cannot be formed at this level. The parent scores in pass B
    # (see _finalize_perm_scores) once every gene's rows are written.
    #
    # Each record is already per side (a "both" position is two records, each with its own
    # `side`), so the detail rows carry `side` directly.
    #
    # ── CT_POSTPROC cluster trains ────────────────────────────────────────────
    # core.postproc.train_flags flags the detected positions of this gene per (base cycle,
    # caap_group). The flag is written per detail row as `clust` (0/1): sub-pass B0 reads it
    # for the dubious-gene test, and B1/B2 skip clust == 1 rows (_is_clustered) when the
    # trains are removed. Nothing is flagged unless postproc_filter is on.
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
            # `r` is already FOP-domain-pooled (or a single-contrast record) and per side;
            # `side` comes from the record and is the direction key downstream.
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
    """caas_row of one null detail row: its asr_path_score (mirror of scoring_compute.R section 2f).

    The observed side takes the same identity, with no rank factor. One definition for all
    the finalize sub-passes.
    """
    return float(row["asr_path_score"])


def _sanitize_gene_shard(gene: str) -> str:
    """Filesystem-safe stem for a per-gene detail shard.

    Gene ids (Ensembl ids, HGNC symbols) are normally clean; this only neutralizes path
    separators so that a stray one cannot escape the shard directory."""
    return gene.replace(os.sep, "__").replace("/", "__").replace("\\", "__").strip() or "_"


def _is_clustered(row: Dict[str, Any]) -> bool:
    """True for a detail row whose position lies in a cluster train (`clust` = 1)."""
    return int(row.get("clust", 0) or 0) == 1


def iter_detail_rows(detail_path: Path):
    """Yield perm_pos_detail rows (dicts) from either layout:

      * a per-gene shard directory  ``perm_pos_detail/<Gene>.tsv.gz``
      * a single concatenated file  ``perm_pos_detail.tsv.gz``

    Shards are read in sorted-filename order and one shard holds exactly one gene, so the
    rows of a gene stay contiguous, which ``_build_cycle_score_pools`` and
    ``_finalize_perm_scores`` pass B2 rely on. A concatenated file must keep that property.
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

    One streaming read of the detail rows. Per (cycle, caap_group, Gene) it accumulates the
    number of distinct detected positions and whether any row lies in a train (`clust`), then
    applies core.postproc.gene_removal, which calibrates within each (cycle, caap_group) pool.
    Returns the (cycle, caap_group, Gene) units to drop from the null pool; empty when `mode`
    is "none".
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
    """Write the (cycle, caap_group, Gene) units dropped by gene removal (TSV), so the removal
    of each labeling can be audited."""
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

    The gene-level null statistic is size_adj_max = F(max)^n (core.scores), the statistic of
    scoring_compute.R section 4a. F must be the ECDF of the pool the positions of the gene were
    drawn from, so every cycle needs its own pool, as the observed side calibrates against its
    own. Otherwise the null would be calibrated against a different reference than the observed
    score and the FCS p.perm comparison would not be valid.

    Directions are kept separate because scoring_compute.R uses direction-matched reference
    pools for the top and bottom scores.

    Pools are stored exactly (sorted array.array('d')), not as a binned histogram: F is raised
    to the power n, so a small ECDF error is amplified n-fold, and the null scores are heavily
    tied, so a bin would absorb many distinct values. A pool is per cycle (the positions of one
    cycle, not of the whole file), so keeping the exact values of all cycles stays affordable.

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
        # agg is keyed (cyc, pos, side) -> {caap_group: caas_row}. core.scores turns each into
        # the score of the position per side, then into one entry per direction: the directional
        # pools take the score of that side, the "all" pool one entry per position (its best
        # side), so a "both" position is not counted twice.
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
    """Pass B: score the detail rows and aggregate them to gene x cycle statistics.

    Follows the observed pipeline of scoring_compute.R:

        null_row_caas   = asr_path_score
        position score  = mean(null_row_caas) over that position's schemes
        gene x cycle    = size_adj_max over the cycle's positions, per direction
                          (CAAS axis; the ASR axis is the 90th percentile of the position scores)

    Runs as two sub-passes over the detail rows: B1 builds the position-score pool of each
    cycle (_build_cycle_score_pools), B2 scores the genes against it. The second read is
    needed because size_adj_max calibrates a gene against the distribution of its own cycle,
    which is known only after every position of that cycle has been seen.

    The detail rows are gene-contiguous (one shard per gene), so a single streaming pass
    aggregates per gene and holds the rows of one gene at a time.

    Outputs in output_dir: gene_cycle_scores.tsv, perm_pos_cycle_caas.tsv.gz,
    perm_pos_sample.tsv and perm_pos_quantiles.tsv.
    """
    import csv as _csv
    import numpy as np

    scores_path = output_dir / "gene_cycle_scores.tsv"
    sample_path = output_dir / "perm_pos_sample.tsv"
    quant_path = output_dir / "perm_pos_quantiles.tsv"
    # Per-(Gene, Position, side, cycle) position score of the null cycle (core.scores), the
    # definition shared by the readers: scoring_compute.R (p.emp, SAM) and the position
    # enrichment. caas_score is empty when no scheme scored the position.
    cycle_caas_path = output_dir / "perm_pos_cycle_caas.tsv.gz"
    cycle_caas_fields = ["Gene", "Position", "side", "cycle", "caas_score", "n_schemes"]

    # Reservoir size per (cycle, scheme). It bounds the sample and the quantile summaries
    # at about K x n_cycles x n_schemes rows regardless of the run size; the full detail
    # rows stay on disk.
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

    # The detail rows always carry a `side` column (a "both" position is two rows, one per
    # side), so positions are keyed by (cyc, pos, side) and the "all" pool keeps the best
    # side per position.

    def _q90(vals) -> float:
        return float(np.percentile(vals, 90)) if vals else 0.0

    def _na(x) -> Any:
        return "NA" if x is None else x

    def _flush(gene: str, pos_scheme: Dict[Any, Dict[str, float]], writer) -> None:
        # pos_scheme is keyed (cyc, pos, side) -> {caap_group: caas_row}. The per-cycle
        # position scores are written first (perm_pos_cycle_caas.tsv.gz).
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
            # The ASR axis is the 90th percentile of the position scores (0 when empty), like
            # the observed gene_caas_asr_* columns; asr_path_score is the caas_row, so it
            # shares the position scores. The CAAS axis is what FCS consumes
            # (caas_corStat_byrank). A cycle with no scored position in a direction has no
            # score there (NA); scoring_caas_perms.R fills those cells with 0 when it builds
            # the dense genes x cycles matrix.
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

            rc = asr  # caas_row is the asr_path_score itself

            skey = (cyc, pos, row.get("side") or "none")
            pos_scheme.setdefault(skey, {})[grp] = rc

            # Reservoir sample stratified by (cycle, scheme): every cycle contributes up to
            # the same K rows however many detections it produced, so each cycle is equally
            # represented and per-cycle summaries can be formed. Sampling per gene instead
            # would weight genes by the number of cycles that hit them.
            key = (cyc, grp)
            seen[key] = seen.get(key, 0) + 1
            res = reservoirs.setdefault(key, [])
            item = (gene, pos, asr)
            if len(res) < sample_per_cycle_group:
                res.append(item)
            else:
                j = rng.randrange(seen[key])
                if j < sample_per_cycle_group:
                    res[j] = item
            n_rows += 1

        if current_gene is not None:
            _flush(current_gene, pos_scheme, writer_scores)

    # ── Sample + quantile summaries ───────────────────────────────────────────
    sample_fields = ["Gene", "Position", "caap_group", "cycle", "asr_path_score"]
    with open(sample_path, "w", newline="") as f_sample:
        writer_sample = _csv.DictWriter(f_sample, fieldnames=sample_fields, delimiter="\t")
        writer_sample.writeheader()
        for (cyc, grp), res in reservoirs.items():
            writer_sample.writerows({
                "Gene": g, "Position": p, "caap_group": grp, "cycle": cyc,
                "asr_path_score": a,
            } for (g, p, a) in res)

    quant_levels = [5, 10, 25, 50, 75, 90, 95]
    quant_fields = (["cycle", "caap_group", "metric", "n_sampled", "n_total", "mean"]
                    + [f"q{q}" for q in quant_levels])
    with open(quant_path, "w", newline="") as f_quant:
        writer_quant = _csv.DictWriter(f_quant, fieldnames=quant_fields, delimiter="\t")
        writer_quant.writeheader()
        for (cyc, grp), res in reservoirs.items():
            if not res:
                continue
            vals = np.asarray([r[2] for r in res], dtype=float)
            rec = {"cycle": cyc, "caap_group": grp, "metric": "asr_path_score",
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
    """Genome-wide CAAS permulation null: load the ASR once per gene, replay N permuted
    labelings, and score them as the observed pipeline scores itself.

    Two passes, because the score of a gene is calibrated against the genome-wide pool of
    position scores of its cycle, while a worker sees a single gene:

      Pass A (parallel, one gene per task, or per chunk of cycles for a large gene) replays
        the labelings and streams the raw per-(gene, cycle, position, scheme) detail to
        perm_pos_detail/<Gene>.tsv.gz (one gz shard per gene, since the reduction of a
        gene yields its complete row list).
      Pass B (single process, streaming) re-reads the shards via iter_detail_rows(),
        builds the pool of position scores of each cycle and derives
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
        outside = [g for g in genes if g not in ensembl_genes]
        if outside:
            logger.warning(
                f"[perms] {len(outside)} of {len(genes)} genes with null hits are not in the Ensembl "
                f"list and are left out of the null, e.g. {sorted(outside)[:5]}"
            )
        genes = [g for g in genes if g in ensembl_genes]

    # FOP mirror: with a fop_pairs_file the resample directory carries "<base>~H<m>"
    # hypothesis labelings. build_cycle_inputs resolves one labeling per hypothesis tag;
    # _perms_worker_replay then domain-pools them back to one record per base cycle. The
    # per-(hypothesis, domain) PSS weights are parsed once here.
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

    # One gz shard per gene rather than one monolithic file: pass A yields one gene's
    # complete row list at a time, and every later re-aggregation (CAAS_CORE_MERGE, the
    # same finalizer) then needs the memory of one gene only and reads shard by shard.
    # `iter_detail_rows()` is the matching reader.
    detail_dir = output_dir / "perm_pos_detail"
    detail_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "perm_pos_detail.manifest.tsv"

    # The detail rows carry `side` (a "both" position is two rows, one per side).
    detail_fields = ["Gene", "cycle", "Position", "caap_group", "asr_path_score",
                     "n_detected", "clust", "side"]

    # CT_POSTPROC filtering of the null candidate pool (cluster trains and dubious/extreme
    # genes). Off by default; the Nextflow layer turns it on (params.caas_perms_postproc)
    # so that the null follows the same filters as the observed filtered_discovery.tsv.
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

    # Optional untrimmed coordinates: with a MAP directory, the cluster trains of each gene
    # measure their span in untrimmed alignment columns (core.columns). A gene without a MAP
    # file keeps the trimmed coordinates.
    train_map_index: Optional[Dict[str, Optional[str]]] = None
    genes_without_map: List[str] = []
    if postproc_filter and train_map_dir:
        from src.core.columns import index_files
        train_map_index = index_files(train_map_dir, train_map_suffix)
        logger.info(f"[perms] cluster trains in untrimmed coordinates: {len(train_map_index)} MAP files in {train_map_dir}")

    # n_cycles_total is the same for every gene (all genes replay the same cycle_tags),
    # so it is computed once here.
    if fop_pairs is not None:
        from src.convergence.fop_pool import base_cycle as _bc_global
        n_cycles_total = len({_bc_global(c) for c in cycle_tags})
    else:
        n_cycles_total = len(cycle_tags)

    if chunk_threshold is None:
        chunk_threshold = int(os.environ.get("CAAS_PERMS_CHUNK_THRESHOLD", "5000"))
    if chunk_target_size is None:
        chunk_target_size = int(os.environ.get("CAAS_PERMS_CHUNK_TARGET_SIZE", "1000"))

    # Chunking: a gene whose size proxy (gene_sizes, the row-count tally that
    # disambiguation_perms_main.py's _genes_from_export computes for the LPT ordering)
    # exceeds chunk_threshold has its replay split into base-cycle-respecting sub-chunks,
    # each dispatched as its own worker task. This bounds the wall time behind a single
    # oversized gene and the memory of a worker, which would otherwise hold the whole
    # multi-thousand-cycle result set of the gene. Below the threshold (or without
    # gene_sizes) a gene gets one chunk with all of cycle_tags. All genes replay the same
    # cycle_tags, so the split points are identical and computed once.
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

    # LPT (longest processing time first): the chunk with the largest estimated workload is
    # dispatched first. `genes` arrives ordered by size (largest first) from
    # _genes_from_export; sorting the flattened chunk list keeps that order when a gene is
    # split into chunks of about equal size.
    tasks.sort(key=lambda t: -t[2])

    args_generator = (
        (gene, alignment_dir, tree_file, taxid_mapping_path, asr_model,
         asr_cache_dir, posterior_threshold,
         chunk, cycle_labelings, perm_discovery_file, None, fop_pairs)
        for gene, chunk, _est in tasks
    )

    # ── Pass A: replay labelings and stream the per-gene detail shards ────────
    cycles_seen: Set[str] = set()
    pool = mp.Pool(
        processes=effective_workers, maxtasksperchild=maxtasks,
        initializer=init_worker, initargs=(1, None),
    )
    n_genes = 0
    n_detail_rows = 0
    manifest_rows: List[Tuple[str, int]] = []
    # Chunk bookkeeping: the pooled results of a gene accumulate here until all its chunks
    # (chunks_per_gene[gene]) have arrived; _perms_worker_finalize then runs once and the
    # entry is dropped. Only the genes with chunks in flight (at most effective_workers) are
    # held partially.
    pending_pooled: Dict[str, List[Tuple[str, List[Any]]]] = {}
    received: Dict[str, int] = {}
    failed_chunks: Dict[str, int] = {}
    try:
        # chunksize=1: each item is the replay of one gene chunk (seconds of work), not
        # the many-cheap-tasks shape chunksize>1 is meant for. A larger chunksize bundles
        # tasks into one unit handed to a single worker, so with few tasks (for example
        # 20 genes and chunksize=10) only as many workers as bundles ever receive work.
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

    # ── Pass B: score each cycle against its own pool, aggregate ──────────────

    # ── Sub-pass B0: cycle-aware dubious/extreme gene removal ─────────────────
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

    # _finalize_perm_scores aggregates the detail shards, which are keyed by base cycle.
    # Under the FOP mirror the `cycle_tags` of build_cycle_inputs are "<base>~H<m>"
    # hypothesis tags, but _perms_worker_replay domain-pools them to one record per base
    # cycle before any detail row is written. The finalizer needs the base-collapsed list:
    # with the hypothesis tags, the `for cyc in cycle_tags` loop of _flush would never match
    # the `cycle` column of the detail rows and every gene x cycle row would be a false
    # zero, leaving gene_cycle_scores.tsv and the FCS p.perm degenerate.
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
