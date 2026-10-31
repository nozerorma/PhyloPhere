#!/usr/bin/env python3
# gene_wrapper.py — Permulation null per gene: replay of the labelings over each gene's ASR, reduced to tile tables.
# PhyloPhere | subworkflows/CT_DISAMBIGUATION/local/src/utils/

"""
The permulation null per gene: replay of the labelings over each gene's ASR, reduced to the tables of a tile.

`process_all_genes_perms` replays the labelings over each gene's alignment in process (src.utils.perm_stream), loads its
ASR once, scores every labeling that has a detection and reduces each cycle as it is produced to a unit row and the
scored positions of the cycle (src.utils.null_tables). The tables of the genes of the tile are written together as
``null_tables/`` (and ``null_tables_strict/``, the strict pooling of the same replay). The gene and position statistics
that need every gene of a cycle (gene removal, the reference pools) are computed later, once per block of cycles, by
src.utils.null_merge. Also the per-gene record converter the master rows are built from.

The null gives CAAS a genome-wide excess null for FCS pathway enrichment: N permuted phenotype labelings are replayed
through the full position pool and scored on the same ASR posteriors as the observed run. The posteriors depend only
on alignment and tree (not on the phenotype), so they are loaded once per gene and reused for every labeling. The
scoring goes through the same code as the observed run, so the null is calibrated against the observed by construction.

Imported by: disambiguation_perms_main.py (process_all_genes_perms), src/core/observed.py (convert_convergence_result_to_dict)
"""

import logging
import multiprocessing as mp
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set, Any

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "src"))

from src.utils.concurrency import plan_concurrency, init_worker, safe_imap_unordered
from src.core.driver import load_gene_context
from src.utils.null_tables import CycleReducer, concat_tables, gene_tables, write_tables
from src.utils.perm_stream import ReplayContext, ReplaySettings, stream_gene

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
        "multi_hypothesis": multi_hypothesis,
        # Hypotheses that drove at least one changed domain on THIS SIDE: the hypothesis-provenance
        # column the scoring reads (scoring_compute.R aggregates it by name). It is never blanked
        # when several hypotheses collide on the position.
        "participating_hypotheses": getattr(result, "participating_hypotheses", None) or "",
        # Number of hypotheses (M) harvested for this position.
        "n_hypotheses": getattr(result, "n_hypotheses", None),
        "n_species": getattr(result, "n_species", None),
        # Cross-hypothesis support tallies for tag and caas, which above keep the
        # value of the first row.
        "tag_support": getattr(result, "tag_support", "") or "",
        "caas_support": getattr(result, "caas_support", "") or "",
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



# ── Permulation null: readers and the per-gene replay ─────────────────────────


def _read_resample_labelings(resample_dir: str) -> Dict[str, Tuple[List[str], List[str]]]:
    """Every resample labeling -> (fg_species, bg_species) by cycle tag (core.labelings.read_labelings)."""
    from src.core.labelings import read_labelings

    return {t: (list(l.fg), list(l.bg)) for t, l in read_labelings(resample_dir).items()}


def _cycle_tags(resample_dir: str, cycles: Optional[List[str]] = None) -> List[str]:
    """The cycle tags to replay: `cycles` when given, else every labeling of `resample_dir` in tag order."""
    return list(cycles) if cycles else sorted(_read_resample_labelings(resample_dir))


# Per-process state of the replay pool, set by _init_replay_worker (a worker started by spawn or forkserver
# receives it through the initializer arguments, once, and parses the labelings of its block itself).
_REPLAY: Dict[str, Any] = {}


def _init_replay_worker(config: Dict[str, Any]) -> None:
    """Pool initializer: thread environment and logging, then the block's labelings and the CAAS kernel."""
    init_worker(1, None)
    _REPLAY.clear()
    _REPLAY.update(config)
    _REPLAY["replay"] = ReplayContext(
        config["settings"], config["resample_dir"], config["cycle_tags"], config["fop_pairs_file"])


def _replay_gene_task(task: Tuple[str, str, Optional[Dict[int, int]]]):
    """Replay of one gene: stream the labelings over its alignment and reduce every cycle to its unit and scores.

    Returns (gene, status, tables, strict tables): status is None when the gene's ASR context could not be loaded (the
    gene is left out of the null, as it is left out of the observed) and 0 for a gene without a detection; otherwise it
    is the number of units. `tables` and `strict tables` are (units, scores) pyarrow tables, None without a detection
    (the strict ones also when the strict pooling is off). A kernel disagreement is not a gene failure: it propagates
    and stops the run.
    """
    gene, alignment_path, train_columns = task
    st = _REPLAY
    replay: ReplayContext = st["replay"]
    try:
        def load_ctx():
            return load_gene_context(
                gene, st["alignment_dir"], st["tree_file"], st["taxid_mapping_path"],
                st["asr_model"], st["asr_cache_dir"], st["posterior_threshold"], None,
                alignment_path=alignment_path, alignment_format=st["settings"].ali_format)

        reducer = CycleReducer(st["postproc_filter"], st["remove_clusters"], st["clust_minlen"], st["clust_maxcaas"], train_columns, st.get("clust_min_positions"),
                               universe=st["universe"][gene])
        # the strict pooling is scored over the strict universe (no conserved pair allowed)
        strict_reducer = CycleReducer(st["postproc_filter"], st["remove_clusters"], st["clust_minlen"], st["clust_maxcaas"], train_columns, st.get("clust_min_positions"),
                                      universe=st["universe_strict"][gene]) if st["strict_output"] else None
        streamed = stream_gene(gene, alignment_path, replay, load_ctx, st["posterior_threshold"], st["strict_output"], reduce=reducer, reduce_strict=strict_reducer)
        if streamed is None:
            return (gene, None, None, None)
        pooled, pooled_strict = streamed
        if not pooled:
            return (gene, 0, None, None)
        tables = gene_tables(gene, [unit for _, unit in pooled])
        strict_tables = gene_tables(gene, [unit for _, unit in pooled_strict]) if pooled_strict else None
        return (gene, len(pooled), tables, strict_tables)
    except replay.pr.KernelDisagreement:
        raise
    except Exception as e:
        logger.error(f"[perms] replay worker failed for {gene}: {e}", exc_info=True)
        return (gene, None, None, None)


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


def process_all_genes_perms(
    genes: List[str],
    gene_alignments: Dict[str, str],
    alignment_dir: str,
    tree_file: str,
    resample_dir: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    workers: Optional[int],
    output_dir: Path,
    replay: ReplaySettings,
    ensembl_genes_file: Optional[str] = None,
    cycles: Optional[List[str]] = None,
    max_tasks_per_child: Optional[int] = None,
    fop_pairs_file: Optional[str] = None,
    clust_minlen: int = 3,
    clust_maxcaas: float = 0.7,
    postproc_filter: bool = False,
    remove_clusters: bool = True,
    train_map_dir: Optional[str] = None,
    train_map_suffix: str = ".map.tsv",
    strict_output: bool = False,
    clust_min_positions: Optional[int] = None,
    universe: Optional[Dict[str, Tuple[int, ...]]] = None,
    universe_strict: Optional[Dict[str, Tuple[int, ...]]] = None,
) -> Path:
    """Replay of the labelings of a tile over each gene, written as the tile tables of the null.

    `genes` are the gene keys to replay, `gene_alignments` the alignment file of each; a gene is replayed against the
    labelings of `resample_dir` (those of `cycles` when given) with the thresholds of `replay`. One task per gene streams
    the labelings over the gene's alignment (src.utils.perm_stream: the CAAS kernel runs in the worker, one range of
    whole base cycles at a time) and reduces each scored cycle to its unit row and scored positions
    (src.utils.null_tables); the parent writes the tables of all the genes together.

    Outputs:
      - output_dir/null_tables/{units,scores}.parquet and meta.json
      - output_dir/null_tables_strict/ (with strict_output): the same replay pooled without the hypotheses whose
        discovery row has conserved pairs; the replay, ASR and scoring are shared with the tables above.
    With `postproc_filter` the cluster trains of each cycle are flagged and, with `remove_clusters`, their rows are left
    out of the scores; the gene-removal thresholds and the reference pools are the merge's.

    `universe` ({gene: positions of the positional universe}, src.core.universe) is required and must cover every gene:
    only its positions are scored, before the gene maximum of each cycle is taken, and its identity per gene is part of
    the tile provenance (meta.json) so a null cannot be reused with another universe. With `strict_output` the strict
    pooling is scored over `universe_strict` (the same test with no conserved pair allowed) and records that one.
    """
    if universe is None or (strict_output and universe_strict is None):
        raise ValueError("[perms] the positional universe is required: the null scores only its positions")
    missing_universe = [g for g in genes if g not in universe or (strict_output and g not in universe_strict)]
    if missing_universe:
        raise ValueError(f"[perms] genes without a positional universe, e.g. {sorted(missing_universe)[:5]}")
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

    # FOP mirror: with a fop_pairs_file the labelings are "<base>~H<m>" hypotheses; the replay domain-pools them back
    # to one record per base cycle with the per-(hypothesis, domain) PSS weights of that file.
    fop_mirror = bool(fop_pairs_file)
    if fop_mirror:
        logger.info(f"[perms] FOP mirror ON: design pairs of {fop_pairs_file}")

    cycle_tags = _cycle_tags(resample_dir, cycles)
    if not cycle_tags:
        raise RuntimeError("[perms] no usable cycles (check the resample dir)")
    # Validate in the parent: initializer exceptions otherwise cause Pool to
    # respawn failing workers indefinitely.
    validated_replay = ReplayContext(replay, str(resample_dir), cycle_tags,
                                     fop_pairs_file if fop_mirror else None)
    del validated_replay
    effective_workers = min(effective_workers, max(1, len(genes)))
    logger.info(
        f"[perms] replaying {len(cycle_tags)} cycles over {len(genes)} genes "
        f"with {effective_workers} workers"
    )

    if max_tasks_per_child is not None:
        maxtasks = int(max_tasks_per_child)
    else:
        maxtasks = int(os.environ.get("CAAS_MAX_TASKS_PER_CHILD", "50"))

    # Optional untrimmed coordinates: with a MAP directory, the cluster trains of each gene
    # measure their span in untrimmed alignment columns (core.columns). A gene without a MAP
    # file keeps the trimmed coordinates.
    train_map_index: Optional[Dict[str, Optional[str]]] = None
    genes_without_map: List[str] = []
    if postproc_filter and train_map_dir:
        from src.core.columns import index_files
        train_map_index = index_files(train_map_dir, train_map_suffix, genes)
        logger.info(f"[perms] cluster trains in untrimmed coordinates: {len(train_map_index)} MAP files in {train_map_dir}")

    # The genes are dispatched in the order given: the caller lists the largest first (longest processing time
    # first), so the biggest gene does not start late and leave workers idle behind it.
    tasks = (
        (gene, str(gene_alignments[gene]), _gene_train_columns(train_map_index, gene, train_map_suffix, genes_without_map))
        for gene in genes
    )
    worker_config = {
        "settings": replay, "resample_dir": str(resample_dir), "cycle_tags": cycle_tags,
        "fop_pairs_file": fop_pairs_file if fop_mirror else None,
        "alignment_dir": alignment_dir, "tree_file": tree_file, "taxid_mapping_path": taxid_mapping_path,
        "asr_model": asr_model, "asr_cache_dir": asr_cache_dir, "posterior_threshold": posterior_threshold,
        "postproc_filter": postproc_filter, "remove_clusters": remove_clusters,
        "clust_min_positions": clust_min_positions, "clust_minlen": clust_minlen, "clust_maxcaas": clust_maxcaas, "strict_output": strict_output,
        "universe": {g: universe[g] for g in genes},
        "universe_strict": {g: universe_strict[g] for g in genes} if strict_output else None,
    }

    cycles_seen: Set[str] = set()
    parts: List[Tuple[Any, Any]] = []
    strict_parts: List[Tuple[Any, Any]] = []
    pool = mp.Pool(
        processes=effective_workers, maxtasksperchild=maxtasks,
        initializer=_init_replay_worker, initargs=(worker_config,),
    )
    n_genes = 0
    n_units = 0
    try:
        # chunksize=1: each item is the replay of one gene (seconds to minutes of work), not the many-cheap-tasks
        # shape chunksize>1 is meant for.
        for _gene, status, tables, strict_tables in safe_imap_unordered(pool, _replay_gene_task, tasks, chunksize=1):
            if status is None:
                raise RuntimeError(f'Null replay failed for {_gene}; incomplete nulls cannot be published')
            if not status:
                continue
            parts.append(tables)
            if strict_tables is not None:
                strict_parts.append(strict_tables)
            cycles_seen.update(tables[0]["cycle"].to_pylist())
            n_units += status
            n_genes += 1
    except BaseException:
        pool.terminate()
        raise
    else:
        pool.close()
    finally:
        pool.join()

    from src.core.scoring_provenance import scientific, digest
    from src.core.universe import digest as universe_digest
    meta = {"by_labeling": not fop_mirror, "strict": False,
            "scientific": scientific(tree_file,asr_model,posterior_threshold,taxid_mapping_path,fop_pairs_file),
            "cluster": dict(enabled=postproc_filter,remove=remove_clusters,minlen=clust_minlen,maxcaas=clust_maxcaas,min_positions=clust_min_positions),
            "alignments": {g:digest(gene_alignments[g]) for g in genes},
            "universe": {g: universe_digest(g, universe[g]) for g in genes},
            "maps": {g: digest(train_map_index.get(g)) if train_map_index else None for g in genes}}

    write_tables(output_dir / "null_tables", *concat_tables(parts), meta=meta)
    if strict_output:
        strict_meta = {**meta, "strict": True,
                       "universe": {g: universe_digest(g, universe_strict[g]) for g in genes}}
        write_tables(output_dir / "null_tables_strict", *concat_tables(strict_parts), meta=strict_meta)
        logger.info(f"[perms] strict pooling: {len(strict_parts)} genes -> null_tables_strict/")

    logger.info(
        f"[perms] replay done: {n_genes} genes, {n_units} (gene, cycle) units "
        f"across {len(cycles_seen)} cycles -> null_tables/"
    )
    if n_genes < len(genes):
        logger.info(
            f"[perms] {n_genes}/{len(genes)} genes with per-cycle CAAS — "
            f"{len(genes) - n_genes} contributed nothing (no CAAS survived any replayed "
            f"labeling, all alignments and ASR contexts were valid)."
        )
    if train_map_index is not None and genes_without_map:
        logger.warning(f"[perms] {len(genes_without_map)} genes have no MAP file and keep trimmed-coordinate "
                       f"trains, e.g. {sorted(genes_without_map)[:5]}")
    if n_units == 0:
        logger.error(
            "[perms] the replay produced ZERO units — the permulation null is empty. "
            "Downstream caas_perms.rds / FCS p.perm will be degenerate. Check the "
            "labelings, the alignments and the ASR cache."
        )
    return output_dir
