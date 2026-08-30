#!/usr/bin/env python3
"""
Gene Processing Wrapper with Multiprocessing

Orchestrates parallel processing of genes using existing single_gene_pipeline
logic. Streams lightweight per-position result dicts into an aggregation SQLite DB,
then exports CSV/JSON from the DB.

Author: ASR Integration
Date: 2025-12-03 (revised 2025-12-09)
"""

import bisect
import gzip
import logging
import multiprocessing as mp
import os
import random
import sys
import threading
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set, Any

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "single_gene_pipeline"))
sys.path.insert(0, str(project_root / "src"))

from src.convergence.disambiguate_single import analyze_gene_disambiguation
from src.asr.asr_single import (
    load_alignment_and_mappings,
    load_and_match_tree,
    run_asr_pipeline,
    load_precomputed_asr,
)

from src.phylo.tree_utils import build_tree_node_mapping, extract_tip_labels
from src.utils.concurrency import plan_concurrency, init_worker, codeml_slot
from src.data.loaders import list_gene_caas_positions, list_gene_caas_entries, normalize_amino_list
from src.utils.io_utils import find_gene_alignment
from src.data.models import CAASPosition

from src.utils.disambiguation_db import (
    init_db,
    get_connection,
    insert_gene_alignment,
    insert_result,
)

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

            if not hasattr(ns, "pvalue") and "pvalue" in result:
                ns.pvalue = result.get("pvalue")

            if not hasattr(ns, "pvalue_boot") and "pvalue_boot" in result:
                ns.pvalue_boot = result.get("pvalue_boot")

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
        "is_significant": bool(getattr(result, "is_significant", False)),
        "pvalue": getattr(result, "pvalue", None),
        "pvalue_boot": getattr(result, "pvalue_boot", None),
        "caap_group": getattr(result, "caap_group", "US"),
        "amino_encoded": getattr(result, "amino_encoded", ""),
        "is_conserved_meta": bool(getattr(result, "is_conserved_meta", False)),
        "conserved_pair": getattr(result, "conserved_pair", ""),
        "sig_hyp": getattr(result, "sig_hyp", None),
        "sig_perm": getattr(result, "sig_perm", None),
        "multi_hypothesis": multi_hypothesis,
    }

    # Node mapping
    node_mapping = getattr(result, "node_mapping", None) or {}
    if isinstance(node_mapping, dict) and node_mapping:
        result_dict["all_mrca_node"] = node_mapping.get("mrca_contrast")

        focal_nodes = []
        focal_nodes_raw = node_mapping.get("focal_nodes")
        if isinstance(focal_nodes_raw, (list, tuple)):
            focal_nodes = list(focal_nodes_raw)
        else:
            # derive from focal_1, focal_2 ...
            candidates = []
            for k, v in node_mapping.items():
                if isinstance(k, str) and k.startswith("focal_"):
                    try:
                        idx = int(k.split("_")[1])
                    except Exception:
                        idx = 0
                    candidates.append((idx, v))
            candidates.sort(key=lambda x: x[0])
            focal_nodes = [v for _, v in candidates]

        for idx, focal_id in enumerate(focal_nodes, 1):
            result_dict[f"mrca_{idx}_node"] = focal_id

    # Node state details
    nsd = getattr(result, "node_state_details", None) or {}
    if isinstance(nsd, dict) and nsd:
        result_dict["all_mrca_state"] = nsd.get("mrca_contrast")
        result_dict["all_mrca_posterior"] = nsd.get("mrca_contrast_prob")

        focal_states = nsd.get("focal_states", []) or []
        focal_probs = nsd.get("focal_probs", []) or []
        for idx in range(1, len(focal_states) + 1):
            state = focal_states[idx - 1] if idx - 1 < len(focal_states) else None
            prob = focal_probs[idx - 1] if idx - 1 < len(focal_probs) else None
            result_dict[f"mrca_{idx}_state"] = state
            result_dict[f"mrca_{idx}_posterior"] = prob

    # Pattern classification
    result_dict["convergence_type"] = getattr(result, "convergence_type", None)

    # Change tracking
    result_dict["change_top"] = getattr(result, "change_top", "no_change")
    result_dict["change_bottom"] = getattr(result, "change_bottom", "no_change")
    result_dict["change_side"] = getattr(result, "change_side", "none")

    # ASR path score (unified ASR/convergence/parallel signal) + per-pair detail
    result_dict["asr_path_score"] = getattr(result, "asr_path_score", None)
    result_dict["independence"] = getattr(result, "independence", None)
    result_dict["mrca_diversity"] = getattr(result, "mrca_diversity", None)
    result_dict["derived_agreement"] = getattr(result, "derived_agreement", None)
    result_dict["conservation_gate"] = getattr(result, "conservation_gate", None)
    result_dict["core"] = getattr(result, "core", None)
    pair_path_scores = getattr(result, "pair_path_scores", None) or {}
    pair_path_contam = getattr(result, "pair_path_contaminated", None) or {}
    if pair_path_scores or pair_path_contam:
        # Fresh ConvergenceResult: per-pair dicts present → flatten to columns.
        for pid, pscore in pair_path_scores.items():
            result_dict[f"mrca_{pid}_path_score"] = pscore
        for pid, contam in pair_path_contam.items():
            result_dict[f"mrca_{pid}_contaminated"] = bool(contam)
    else:
        # Round-trip case (reloaded from aggregation DB): the per-pair scores
        # were already flattened onto the input as mrca_N_path_score /
        # mrca_N_contaminated attributes. Carry those flat keys through.
        src = vars(result) if hasattr(result, "__dict__") else {}
        for k, v in src.items():
            if k.startswith("mrca_") and (
                k.endswith("_path_score") or k.endswith("_contaminated")
            ):
                result_dict[k] = v

    # Keep rich payloads for downstream viz/reporting
    if getattr(result, "pair_details", None):
        result_dict["pair_details"] = getattr(result, "pair_details")
    if getattr(result, "node_mapping", None):
        result_dict["node_mapping"] = getattr(result, "node_mapping")
    if getattr(result, "node_state_details", None):
        result_dict["node_state_details"] = getattr(result, "node_state_details")

    return result_dict


def merge_multi_hypothesis_results(
    results_group: List,
    alignment=None,
    seq_by_id: Optional[Dict] = None,
    seq_by_species: Optional[Dict] = None,
    trait_pairs: Optional[Dict[int, List[Tuple[str, str]]]] = None,
    taxid_to_species: Optional[Dict] = None,
) -> Dict:
    """
    Merge multiple ConvergenceResult rows for same (gene, msa_pos) into one dict.

    This function is intentionally side-effect free: it does not mutate
    the input results.
    """
    base = results_group[0]
    tags = [
        str(getattr(r, "tag", "")) for r in results_group if getattr(r, "tag", None)
    ]
    multi_hypothesis = ",".join(sorted(set(tags))) if tags else None

    merged = convert_convergence_result_to_dict(
        base,
        multi_hypothesis=multi_hypothesis,
        alignment=alignment,
        seq_by_id=seq_by_id,
        seq_by_species=seq_by_species,
        trait_pairs=trait_pairs,
        taxid_to_species=taxid_to_species,
    )

    # Merge CAAS labels
    caass = [
        getattr(r, "caas", None) for r in results_group if getattr(r, "caas", None)
    ]
    if caass:
        merged["caas_merged"] = "; ".join(sorted(set(str(x) for x in caass)))

    # Merge boolean/meta flags across hypothesis rows
    merged["is_conserved_meta"] = any(
        bool(getattr(r, "is_conserved_meta", False)) for r in results_group
    )
    path_scores = [
        getattr(r, "asr_path_score", None)
        for r in results_group
        if getattr(r, "asr_path_score", None) is not None
    ]
    if path_scores:
        merged["asr_path_score"] = max(path_scores)

    groups = [
        str(getattr(r, "caap_group", "")).strip()
        for r in results_group
        if getattr(r, "caap_group", None)
    ]
    if groups:
        merged["caap_group"] = ",".join(sorted(set(g for g in groups if g)))

    encoded = [
        str(getattr(r, "amino_encoded", "")).strip()
        for r in results_group
        if getattr(r, "amino_encoded", None)
    ]
    if encoded:
        merged["amino_encoded"] = ",".join(sorted(set(e for e in encoded if e)))

    return merged


def process_single_gene(
    gene: str,
    alignment_dir: str,
    tree_file: str,
    caas_metadata_path: str,
    trait_file_path: str,
    taxid_mapping_path: Optional[str],
    asr_mode: str,
    asr_model: str,
    asr_cache_dir: Optional[str],
    posterior_threshold: float,
    convergence_mode: str,
    threads_per_gene: int,
    run_diagnostics: bool,
    output_dir: Path,
    db_queue: Optional[Any] = None,
    ensembl_genes: Optional[Set[str]] = None,
) -> Tuple[str, Optional[Path]]:

    try:
        alignment_path = find_gene_alignment(Path(alignment_dir), gene, ensembl_genes)
        if not alignment_path:
            logger.warning(f"No alignment found for {gene}, skipping")
            return (gene, None)

        caas_positions = list_gene_caas_positions(Path(caas_metadata_path), gene)
        if not caas_positions:
            logger.debug(f"No CAAS positions for {gene}")
            return (gene, None)

        logger.debug(f"Processing gene: {gene} ({len(caas_positions)} CAAS positions)")

        alignment_data = load_alignment_and_mappings(
            alignment_path,
            Path(taxid_mapping_path) if taxid_mapping_path else None,
            gene_name=gene,
        )

        tree_data = load_and_match_tree(
            Path(tree_file),
            alignment_data,
            Path(taxid_mapping_path) if taxid_mapping_path else None,
        )

        node_posteriors = None
        rst_file = None
        paml_tree_file = None

        if asr_mode == "precomputed" and asr_cache_dir:
            from src.asr.asr_single import SingleGeneASRConfig

            asr_config = SingleGeneASRConfig(
                alignment_path=alignment_path,
                tree_path=Path(tree_file),
                model=asr_model,
                posterior_threshold=posterior_threshold,
                output_dir=Path(asr_cache_dir),
            )
            try:
                node_posteriors = load_precomputed_asr(gene, asr_config, alignment_data)
                if node_posteriors:
                    rst_file = getattr(node_posteriors, "rst_file", None)
                    paml_tree_file = getattr(node_posteriors, "tree_file", None)
            except FileNotFoundError as e:
                logger.debug(f"Precomputed ASR not found for {gene}: {e}")

        elif asr_mode == "compute":
            gene_output_dir = Path(asr_cache_dir) if asr_cache_dir else output_dir / "asr"
            gene_output_dir.mkdir(parents=True, exist_ok=True)

            from src.asr.asr_single import SingleGeneASRConfig

            asr_config = SingleGeneASRConfig(
                alignment_path=alignment_path,
                tree_path=Path(tree_file),
                model=asr_model,
                posterior_threshold=posterior_threshold,
                output_dir=gene_output_dir,
                threads=threads_per_gene,
            )

            with codeml_slot():
                node_posteriors = run_asr_pipeline(
                    gene,
                    asr_config,
                    skip_if_exists=True,
                    alignment_data=alignment_data,
                    tree_data=tree_data,
                )
            if node_posteriors:
                rst_file = getattr(node_posteriors, "rst_file", None)
                paml_tree_file = getattr(node_posteriors, "tree_file", None)

        # Rebuild tree_data with PAML-labeled tree to align node IDs
        if rst_file and paml_tree_file and Path(paml_tree_file).exists():
            try:
                ordered_nodes, id_mapping = build_tree_node_mapping(
                    tree_file=Path(paml_tree_file), rst_file=Path(rst_file)
                )
                tree_data.nodes = ordered_nodes
                tree_data.root = ordered_nodes[-1]
                tree_data.node_mapping = id_mapping

                def _tip_taxid(label: str) -> str:
                    return label.split("_")[-1] if "_" in label else label

                tree_data.tip_set = {
                    _tip_taxid(lbl) for lbl in extract_tip_labels(tree_data.root)
                }
                logger.debug(
                    "Updated tree_data from PAML tree for node/posterior alignment"
                )
            except Exception as e:
                logger.warning(f"Could not rebuild tree_data from PAML tree: {e}")

        diag_root = output_dir / "diagnostics" if run_diagnostics else None
        posterior_dump_jsonl = None

        if (
            run_diagnostics
            and node_posteriors
            and getattr(node_posteriors, "posteriors_node", None)
        ):
            try:
                from src.asr.posterior import export_posteriors_to_jsonl

                effective_diag_root = (
                    diag_root if diag_root is not None else (output_dir / "diagnostics")
                )
                posterior_dump_dir = effective_diag_root / "node_dumps"
                posterior_dump_dir.mkdir(parents=True, exist_ok=True)
                posterior_dump_jsonl = (
                    posterior_dump_dir / f"{gene.lower()}_posteriors.jsonl"
                )

                if not posterior_dump_jsonl.exists():
                    posteriors_node = (
                        getattr(node_posteriors, "posteriors_node", {}) or {}
                    )
                    posteriors_mapping = {int(k): v for k, v in posteriors_node.items()}
                    export_posteriors_to_jsonl(
                        posteriors_mapping,
                        posterior_dump_jsonl,
                        (
                            node_posteriors.node_id_map
                            if hasattr(node_posteriors, "node_id_map")
                            else None
                        ),
                    )
                    logger.debug(f"Wrote posterior JSONL to {posterior_dump_jsonl}")
            except Exception as e:
                logger.warning(f"Failed to export posterior JSONL for {gene}: {e}")

        # Filter posteriors to CAAS positions
        filtered_posteriors = None
        if node_posteriors and getattr(node_posteriors, "posteriors_node", None):
            positions_1_based = {int(p) + 1 for p in caas_positions}
            try:
                from src.asr.posterior import parse_paml_rst_node_level

                if getattr(node_posteriors, "rst_file", None) and getattr(
                    node_posteriors, "tree_file", None
                ):
                    parsed = parse_paml_rst_node_level(
                        Path(str(node_posteriors.rst_file)),
                        Path(str(node_posteriors.tree_file)),
                        threshold=posterior_threshold,
                        positions=positions_1_based,
                    )
                    filtered_posteriors = (
                        parsed[0] if isinstance(parsed, tuple) else parsed
                    )
            except Exception:
                filtered_posteriors = None

            if filtered_posteriors is None:
                try:
                    posteriors_node = (
                        getattr(node_posteriors, "posteriors_node", {}) or {}
                    )
                    filtered_posteriors = {
                        int(node_id): {
                            int(pos): probs
                            for pos, probs in site_map.items()
                            if int(pos) in positions_1_based
                        }
                        for node_id, site_map in posteriors_node.items()
                    }
                except Exception:
                    filtered_posteriors = None

        biochem_results, diagnostics = analyze_gene_disambiguation(
            gene=gene,
            alignment_data=alignment_data,
            tree_data=tree_data,
            caas_positions=caas_positions,
            caas_metadata_path=Path(caas_metadata_path),
            trait_file_path=Path(trait_file_path),
            taxid_mapping=alignment_data.species_to_taxid,
            posterior_data=(
                filtered_posteriors
                if filtered_posteriors
                else (node_posteriors.posteriors_node if node_posteriors else None)
            ),
            posterior_threshold=posterior_threshold,
            diagnostics_dir=diag_root,
            convergence_mode=convergence_mode,
            asr_mode=asr_mode,
        )

        if not biochem_results:
            return (gene, None)

        # Prepare JSON-serializable alignment record
        seq_by_id_raw = getattr(alignment_data, "seq_by_id", {}) or {}
        seq_by_species_raw = getattr(alignment_data, "seq_by_species", {}) or {}

        def _to_seq_str(val):
            try:
                if hasattr(val, "seq"):
                    return str(val.seq)
                if isinstance(val, str):
                    return val
                return str(val)
            except Exception:
                return None

        seq_by_id_serial = {k: _to_seq_str(v) for k, v in seq_by_id_raw.items()}
        seq_by_species_serial = {
            k: _to_seq_str(v) for k, v in seq_by_species_raw.items()
        }

        alignment_record = {
            "alignment_path": str(alignment_path) if alignment_path else None,
            "seq_by_id": seq_by_id_serial,
            "seq_by_species": seq_by_species_serial,
            "taxid_to_species": getattr(alignment_data, "taxid_to_species", None),
            "species_to_taxid": getattr(alignment_data, "species_to_taxid", None),
            "alignment_extras": {
                "paml_tree_file": paml_tree_file,
                "posterior_dump_jsonl": (
                    str(posterior_dump_jsonl) if posterior_dump_jsonl else None
                ),
                "node_id_map": (
                    node_posteriors.node_id_map
                    if node_posteriors and hasattr(node_posteriors, "node_id_map")
                    else None
                ),
            },
        }

        if db_queue is not None:
            try:
                db_queue.put(
                    {"type": "gene", "gene": gene, "alignment": alignment_record}
                )

                for r in biochem_results:
                    msa_pos = getattr(r, "position", None)
                    caas_dict = convert_convergence_result_to_dict(
                        r,
                        multi_hypothesis=None,
                        alignment=None,
                        seq_by_id=seq_by_id_serial,
                        seq_by_species=seq_by_species_serial,
                        trait_pairs=None,
                        taxid_to_species=getattr(
                            alignment_data, "taxid_to_species", None
                        ),
                    )

                    db_queue.put(
                        {
                            "type": "result",
                            "gene": gene,
                            "msa_pos": msa_pos,
                            "position": caas_dict.get("position", None),
                            "result": caas_dict,
                        }
                    )

            except Exception as e:
                logger.warning(f"Failed to enqueue results for {gene}: {e}")

        logger.debug(f"  {gene}: {len(biochem_results)} CAAS processed (streamed to DB)")

        # keep worker RAM low
        try:
            import gc

            for name in (
                "alignment_data",
                "tree_data",
                "biochem_results",
                "node_posteriors",
                "seq_by_id_raw",
                "seq_by_species_raw",
                "seq_by_id_serial",
                "seq_by_species_serial",
            ):
                if name in locals():
                    try:
                        del locals()[name]
                    except Exception:
                        pass
            gc.collect()
        except Exception:
            pass

        return (gene, None)

    except FileNotFoundError:
        return (gene, None)

    except Exception as e:
        logger.error(f"Failed to process {gene}: {e}", exc_info=True)
        return (gene, None)


def process_all_genes(
    genes: List[str],
    alignment_dir: str,
    tree_file: str,
    caas_metadata_path: str,
    trait_file_path: str,
    taxid_mapping_path: Optional[str],
    asr_mode: str,
    asr_model: str,
    asr_cache_dir: Optional[str],
    posterior_threshold: float,
    convergence_mode: str,
    threads_per_gene: int,
    workers: Optional[int],
    run_diagnostics: bool,
    output_dir: Path,
    ensembl_genes_file: Optional[str] = None,
    max_tasks_per_child: Optional[int] = None,
    max_codeml: Optional[int] = None,
) -> Tuple[List[Dict], Optional[Dict]]:

    effective_workers, threads_per_gene = plan_concurrency(
        workers, threads_per_gene, logger
    )
    logger.info(
        f"Processing {len(genes)} genes with {effective_workers} workers and {threads_per_gene} threads/gene"
    )

    db_queue = None
    writer_thread = None
    db_path = output_dir / "aggregation.sqlite3"

    # Load Ensembl genes file if provided (no gate)
    from src.data.loaders import load_ensembl_genes

    ensembl_genes: Optional[Set[str]] = None
    try:
        if ensembl_genes_file:
            ensembl_genes = load_ensembl_genes(Path(ensembl_genes_file)) or set()
            logger.info(f"Loaded {len(ensembl_genes)} genes from {ensembl_genes_file}")
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning(f"Failed to load Ensembl genes from {ensembl_genes_file}: {exc}")

    # Optional gate to limit concurrent codeml runs
    codeml_sem = None
    if max_codeml is not None:
        manager = mp.Manager()
        codeml_sem = manager.Semaphore(max(1, int(max_codeml)))

    # Init DB + queue + writer
    init_db(db_path)
    manager = mp.Manager()
    db_queue = manager.Queue()

    def _db_writer(db_path_local, queue):
        conn = get_connection(db_path_local)
        genes_seen = set()
        insert_count = 0
        try:
            while True:
                item = queue.get()
                if item is None:
                    break
                try:
                    if item.get("type") == "gene":
                        gene_name = item.get("gene")
                        if gene_name and gene_name not in genes_seen:
                            insert_gene_alignment(
                                conn, gene_name, item.get("alignment") or {}
                            )
                            genes_seen.add(gene_name)

                    elif item.get("type") == "result":
                        insert_result(
                            conn,
                            item.get("gene"),
                            item.get("msa_pos"),
                            item.get("position"),
                            item.get("result") or {},
                        )
                        insert_count += 1
                        if insert_count % 100 == 0:
                            conn.commit()

                except Exception:
                    conn.rollback()
                    raise
        finally:
            conn.commit()
            conn.close()

    writer_thread = threading.Thread(
        target=_db_writer, args=(db_path, db_queue), daemon=True
    )
    writer_thread.start()

    if max_tasks_per_child is not None:
        maxtasks = int(max_tasks_per_child)
    else:
        maxtasks = int(os.environ.get("CAAS_MAX_TASKS_PER_CHILD", "50"))

    pool = mp.Pool(
        processes=effective_workers,
        initializer=init_worker,
        initargs=(threads_per_gene, codeml_sem),
        maxtasksperchild=maxtasks,
    )

    async_results = []
    try:
        for gene in genes:
            async_results.append(
                pool.apply_async(
                    process_single_gene,
                    (
                        gene,
                        alignment_dir,
                        tree_file,
                        caas_metadata_path,
                        trait_file_path,
                        taxid_mapping_path,
                        asr_mode,
                        asr_model,
                        asr_cache_dir,
                        posterior_threshold,
                        convergence_mode,
                        threads_per_gene,
                        run_diagnostics,
                        output_dir,
                        db_queue,
                        ensembl_genes,
                    ),
                )
            )

        # Force retrieval to surface errors early
        for idx, async_res in enumerate(async_results):
            gene = genes[idx]
            try:
                async_res.get()
            except Exception as e:
                logger.error(f"Gene {gene} failed: {e}", exc_info=True)

    finally:
        pool.close()
        pool.join()

    # Finish DB writer
    db_queue.put(None)
    writer_thread.join()

    # Determine processed genes
    conn = get_connection(db_path)
    try:
        cur = conn.cursor()
        cur.execute("SELECT DISTINCT gene FROM results")
        processed_genes = {row[0] for row in cur.fetchall()}
    finally:
        conn.close()

    skipped_genes = [g for g in genes if g not in processed_genes]
    if skipped_genes:
        diag_dir = output_dir / "diagnostics"
        diag_dir.mkdir(parents=True, exist_ok=True)
        skipped_file = diag_dir / "skipped_genes.txt"
        with open(skipped_file, "w") as f:
            f.write(
                f"# Skipped {len(skipped_genes)} genes (missing alignment or no CAAS positions)\n"
            )
            for g in sorted(skipped_genes):
                f.write(f"{g}\n")
        logger.info(f"Skipped {len(skipped_genes)} genes (see {skipped_file})")

    logger.info(
        f"Results written to DB at {db_path}; processed {len(processed_genes)} genes"
    )

    from src.reporting.disambiguation_writers import export_from_db

    caas_files, summary_json = export_from_db(db_path, output_dir)

    export_info = {
        "db_path": str(db_path),
        "caas_files": [str(p) for p in caas_files],
        "summary_json": str(summary_json),
    }

    return [], export_info


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
    """Index every resample cycle → (fg_species, bg_species) by its cycle tag.

    resample_*.tab is 3-col, no header: cycle_tag, fg_csv, bg_csv (permulations.R).
    A file may hold many cycles (chunk_size rows). Pairing is by index: fg[k] ↔ bg[k]
    is pair k+1 — exactly how the original trait file encodes pairs.
    """
    import csv

    labelings: Dict[str, Tuple[List[str], List[str]]] = {}
    for tab in sorted(Path(resample_dir).glob("resample_*.tab")):
        try:
            with open(tab, "r") as f:
                for row in csv.reader(f, delimiter="\t"):
                    if len(row) < 3:
                        continue
                    tag = row[0].strip()
                    fg = [s for s in row[1].split(",") if s.strip()]
                    bg = [s for s in row[2].split(",") if s.strip()]
                    if tag and fg and bg:
                        labelings[tag] = (fg, bg)
        except Exception as e:
            logger.warning(f"[perms] failed to read {tab}: {e}")
    return labelings


def _write_cycle_trait_file(fg: List[str], bg: List[str], out_path: Path) -> None:
    """Write a cycle's labeling as a 3-col trait file (species, trait, pair).

    trait=1 (high/top) for FG, trait=0 (low/bottom) for BG; pair index = position in
    the resample lists. Consumed verbatim by parse_trait_pairs in disambiguation.
    """
    with open(out_path, "w") as f:
        for k, sp in enumerate(fg, start=1):
            f.write(f"{sp}\t1\t{k}\n")
        for k, sp in enumerate(bg, start=1):
            f.write(f"{sp}\t0\t{k}\n")


def build_cycle_inputs(
    perm_discovery_path: str,
    resample_dir: str,
    work_dir: Path,
    cycles: Optional[List[str]] = None,
) -> Tuple[List[str], Dict[str, str], Dict[str, str]]:
    """Materialize per-cycle trait + caas-metadata files from the resample labelings.
    Skips writing metadata split files to disk since workers read gene files directly.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    trait_dir = work_dir / "cycle_traits"
    trait_dir.mkdir(exist_ok=True)

    labelings = _read_resample_labelings(resample_dir)
    cycle_trait_files: Dict[str, str] = {}
    cycle_meta_files: Dict[str, str] = {}

    target_cycles = cycles if cycles else sorted(labelings.keys())

    for c in target_cycles:
        if c in labelings:
            fg, bg = labelings[c]
            trait_path = trait_dir / f"trait_{c}.tab"
            if not trait_path.exists():
                _write_cycle_trait_file(fg, bg, trait_path)
            cycle_trait_files[c] = str(trait_path)
            cycle_meta_files[c] = ""

    logger.info(f"[perms] prepared trait inputs for {len(target_cycles)} cycles")
    return target_cycles, cycle_trait_files, cycle_meta_files




def _load_gene_asr_context(
    gene: str,
    alignment_dir: str,
    tree_file: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    ensembl_genes: Optional[Set[str]] = None,
) -> Optional[Dict[str, Any]]:
    """Load alignment + tree + precomputed ASR posteriors for one gene ONCE.

    Mirrors the precomputed-ASR load path of process_single_gene (the
    phenotype-invariant part), including the PAML tree rebuild for node alignment.
    Returns a context dict, or None if alignment/ASR unavailable.
    """
    alignment_path = find_gene_alignment(Path(alignment_dir), gene, ensembl_genes)
    if not alignment_path:
        return None

    alignment_data = load_alignment_and_mappings(
        alignment_path,
        Path(taxid_mapping_path) if taxid_mapping_path else None,
        gene_name=gene,
    )
    tree_data = load_and_match_tree(
        Path(tree_file), alignment_data,
        Path(taxid_mapping_path) if taxid_mapping_path else None,
    )

    from src.asr.asr_single import SingleGeneASRConfig

    asr_config = SingleGeneASRConfig(
        alignment_path=alignment_path,
        tree_path=Path(tree_file),
        model=asr_model,
        posterior_threshold=posterior_threshold,
        output_dir=Path(asr_cache_dir),
    )
    try:
        node_posteriors = load_precomputed_asr(gene, asr_config, alignment_data)
    except FileNotFoundError:
        return None
    if not node_posteriors:
        return None

    rst_file = getattr(node_posteriors, "rst_file", None)
    paml_tree_file = getattr(node_posteriors, "tree_file", None)
    if rst_file and paml_tree_file and Path(paml_tree_file).exists():
        try:
            ordered_nodes, id_mapping = build_tree_node_mapping(
                tree_file=Path(paml_tree_file), rst_file=Path(rst_file)
            )
            tree_data.nodes = ordered_nodes
            tree_data.root = ordered_nodes[-1]
            tree_data.node_mapping = id_mapping

            def _tip_taxid(label: str) -> str:
                return label.split("_")[-1] if "_" in label else label

            tree_data.tip_set = {
                _tip_taxid(lbl) for lbl in extract_tip_labels(tree_data.root)
            }
        except Exception as e:
            logger.warning(f"[perms] {gene}: could not rebuild tree_data from PAML tree: {e}")

    return {
        "alignment_data": alignment_data,
        "tree_data": tree_data,
        "node_posteriors": node_posteriors,
    }


# Convergence states that count as "this position changed on this side". Shared by
# the perms worker and the pass-B finalizer so change_side can never drift between
# the two halves of the null computation.
_CHANGE_STATES = ("convergent", "codivergent", "divergent")

# No per-scheme weight any more. scoring_compute.R section 2g aggregates a
# position's schemes with a MEAN of caas_row, not a 0.2-weighted sum, because the
# number of detecting schemes is a biochemical-distance property of the
# substitution rather than evidence strength. The null mirrors that exactly.


def build_percent_rank_lookup(hist_by_cycle: Dict[str, Dict[int, int]]) -> Dict[str, Dict[int, float]]:
    """Exact dplyr::percent_rank over each cycle's genome-wide candidate pool.

    The null's phenotype axis is a rank, matching how the observed side builds its
    own (scoring_compute.R: `1 - percent_rank(pvalue_boot)` over the whole scored
    pool). The rank is what makes the axis usable: the underlying recovery
    statistic is heavily concentrated near zero (mean ~0.03), so `1 - raw` would be
    a near-constant ~0.97, whereas `1 - percent_rank` is uniform on [0, 1] with
    mean 0.50 — the same scale the observed axis occupies, which is the condition
    for the two being comparable at all.

    Each cycle is ranked within its OWN candidate pool, so the pool is the set of
    (gene, position, scheme) triples that cycle discovered. This is what gives the
    axis per-cycle variation: a position's raw LOO value is fixed, but its rank
    depends on the composition of the cycle it appears in.

    `null_pvalue_boot` is monotone in n_detected, so ranking n_detected ranks the
    p-value identically. That collapses the pool to a per-cycle histogram
    (n_detected -> count), a few MB rather than tens of millions of rows, which is
    what allows the rank to be formed without holding the full detail in memory.

    dplyr defines percent_rank(x) = (rank(x, ties="min") - 1) / (n - 1), i.e. the
    fraction of the pool STRICTLY BELOW x. A cumulative sweep over the histogram
    reproduces that exactly, ties included.

    Returns {cycle: {n_detected: percent_rank}}.
    """
    lookup: Dict[str, Dict[int, float]] = {}
    for cycle, hist in hist_by_cycle.items():
        total = sum(hist.values())
        # dplyr yields NaN for a single-element vector; a lone candidate has
        # nothing to be ranked against, so treat it as the bottom of its pool
        # (percent_rank 0 -> phen_score 1), consistent with the strictly-below
        # definition.
        if total <= 1:
            lookup[cycle] = {d: 0.0 for d in hist}
            continue
        denom = float(total - 1)
        cum_below = 0
        per_cycle: Dict[int, float] = {}
        for d in sorted(hist):
            per_cycle[d] = cum_below / denom
            cum_below += hist[d]
        lookup[cycle] = per_cycle
    return lookup


def _perms_worker(
    gene: str,
    alignment_dir: str,
    tree_file: str,
    taxid_mapping_path: Optional[str],
    asr_model: str,
    asr_cache_dir: str,
    posterior_threshold: float,
    convergence_mode: str,
    cycle_tags: List[str],
    cycle_trait_files: Dict[str, str],
    perm_discovery_file: str,
    ensembl_genes: Optional[Set[str]] = None,
) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Replay N labelings over one gene's cached ASR. One gene = one worker task.
    Performs on-the-fly aggregation inside the worker, avoiding outputting huge files.
    """
    try:
        ctx = _load_gene_asr_context(
            gene, alignment_dir, tree_file, taxid_mapping_path,
            asr_model, asr_cache_dir, posterior_threshold, ensembl_genes,
        )
        if ctx is None:
            # No cached ASR for this gene (asr-cache-dir is precomputed-only, no
            # compute-on-miss). Genes that only ever appear under a NULL labeling
            # are never processed by the observed-run disambiguation and so have
            # no cache — they are excluded from the null. Warn rather than drop
            # silently.
            logger.warning(
                f"[perms] no cached ASR for {gene} in {asr_cache_dir} — "
                f"excluded from the permulation null"
            )
            return (gene, [], [])

        alignment_data = ctx["alignment_data"]
        tree_data = ctx["tree_data"]
        node_posteriors = ctx["node_posteriors"]
        full_posteriors = getattr(node_posteriors, "posteriors_node", None)

        axes_only = os.environ.get("CAAS_PERMS_AXES_ONLY", "1") not in ("0", "false", "False")
        per_site_dist_cache: Dict[int, Any] = {} if axes_only else None

        # Load the gene's bootstrap discovery output in memory once
        cycle_to_entries = {}
        gene_file = None
        disc_path = Path(perm_discovery_file)
        if disc_path.is_file():
            with open(disc_path, "r") as f:
                header = f.readline()
                if header:
                    cols = header.rstrip("\n").split("\t")
                    col_indices = {col: idx for idx, col in enumerate(cols)}
                    if "cycle" in col_indices and "gene" in col_indices and "position" in col_indices:
                        cyc_idx = col_indices["cycle"]
                        gene_idx = col_indices["gene"]
                        pos_idx = col_indices["position"]
                        caas_idx = col_indices.get("caas")
                        ae_idx = col_indices.get("amino_encoded")
                        icm_idx = col_indices.get("is_conserved_meta")
                        cp_idx = col_indices.get("conserved_pair")
                        pval_idx = col_indices.get("pvalue")
                        
                        for line in f:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) <= max(cyc_idx, gene_idx, pos_idx):
                                continue
                            if parts[gene_idx].strip() != gene:
                                continue
                            cyc = parts[cyc_idx].strip()
                            if not cyc or cyc not in cycle_tags:
                                continue
                            
                            pos_raw = parts[pos_idx].strip()
                            try:
                                pos0 = int(pos_raw)
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
                                pvalue=float(parts[pval_idx]) if pval_idx is not None and pval_idx < len(parts) and parts[pval_idx] not in ("", "NA") else None,
                                pvalue_boot=None,
                                is_significant=False,
                                caap_group=parts[col_indices["caap_group"]] if "caap_group" in col_indices and col_indices["caap_group"] < len(parts) else "US",
                                amino_encoded=parts[ae_idx] if ae_idx is not None and ae_idx < len(parts) else "",
                                is_conserved_meta=parts[icm_idx] in ("TRUE", "True", "1") if icm_idx is not None and icm_idx < len(parts) else False,
                                conserved_pair=cp,
                            )
                            if cyc not in cycle_to_entries:
                                cycle_to_entries[cyc] = []
                            cycle_to_entries[cyc].append(entry)
        else:
            for p in disc_path.iterdir():
                if p.is_file() and p.name.split(".", 1)[0] == gene:
                    gene_file = p
                    break
            if gene_file and gene_file.exists():
                with open(gene_file, "r") as f:
                    header = f.readline()
                    if header:
                        cols = header.rstrip("\n").split("\t")
                        col_indices = {col: idx for idx, col in enumerate(cols)}
                        if "cycle" in col_indices and "position" in col_indices:
                            cyc_idx = col_indices["cycle"]
                            pos_idx = col_indices["position"]
                            caas_idx = col_indices.get("caas")
                            ae_idx = col_indices.get("amino_encoded")
                            icm_idx = col_indices.get("is_conserved_meta")
                            cp_idx = col_indices.get("conserved_pair")
                            pval_idx = col_indices.get("pvalue")
                            
                            for line in f:
                                parts = line.rstrip("\n").split("\t")
                                if len(parts) <= max(cyc_idx, pos_idx):
                                    continue
                                cyc = parts[cyc_idx].strip()
                                if not cyc or cyc not in cycle_tags:
                                    continue
                                
                                pos_raw = parts[pos_idx].strip()
                                try:
                                    pos0 = int(pos_raw)
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
                                    pvalue=float(parts[pval_idx]) if pval_idx is not None and pval_idx < len(parts) and parts[pval_idx] not in ("", "NA") else None,
                                    pvalue_boot=None,
                                    is_significant=False,
                                    caap_group=parts[col_indices["caap_group"]] if "caap_group" in col_indices and col_indices["caap_group"] < len(parts) else "US",
                                    amino_encoded=parts[ae_idx] if ae_idx is not None and ae_idx < len(parts) else "",
                                    is_conserved_meta=parts[icm_idx] in ("TRUE", "True", "1") if icm_idx is not None and icm_idx < len(parts) else False,
                                    conserved_pair=cp,
                                )
                                if cyc not in cycle_to_entries:
                                    cycle_to_entries[cyc] = []
                                cycle_to_entries[cyc].append(entry)

        all_cycle_results = []
        for cyc in cycle_tags:
            trait_path = cycle_trait_files.get(cyc)
            if not trait_path:
                continue
            caas_entries = cycle_to_entries.get(cyc, [])
            if not caas_entries:
                continue
            try:
                biochem_results, _ = analyze_gene_disambiguation(
                    gene=gene,
                    alignment_data=alignment_data,
                    tree_data=tree_data,
                    caas_positions=[],
                    caas_entries=caas_entries,
                    caas_metadata_path=Path("dummy_path"),
                    trait_file_path=Path(trait_path),
                    taxid_mapping=alignment_data.species_to_taxid,
                    posterior_data=full_posteriors,
                    posterior_threshold=posterior_threshold,
                    diagnostics_dir=None,
                    convergence_mode=convergence_mode,
                    asr_mode="precomputed",
                    axes_only=axes_only,
                    per_site_dist_cache=per_site_dist_cache,
                )
                if biochem_results:
                    all_cycle_results.append((cyc, biochem_results))
            except Exception as e:
                logger.debug(f"[perms] {gene} cycle {cyc} failed: {e}")
                continue

        if not all_cycle_results:
            return (gene, [], [])

        # ── 1. Leave-one-out null_pvalue_boot per (position, scheme) ───────────
        # This is the null-side analogue of the observed pvalue_boot, which
        # modules/boot.py computes as count/cycles: "of the resampling cycles,
        # what fraction still recovered THIS position". Here the resampling axis
        # is phenotype relabeling rather than species bootstrap, but the
        # statistic is the same shape and its empirical distribution matches
        # closely (observed mean 0.048 vs null 0.030 on cancer_complete_bm).
        #
        # LEAVE-ONE-OUT: when scoring cycle i, cycle i must not count toward its
        # own replication evidence, or every cycle's score is contaminated by the
        # very draw it is meant to be an independent sample of. Since a position
        # is only ever scored in the cycles that detected it, k_{-i} is simply
        # n_detected - 1 for every such cycle, so the LOO value is one number per
        # (position, scheme) rather than per (cycle, position, scheme).
        #
        # NOTE the LOO correction is numerically inert once the caller applies
        # percent_rank downstream ((d-1)/(N-1) and d/N are both monotone in d, so
        # they rank identically). It is kept because perm_pos_pval.tsv is read
        # directly as an observed-vs-null crosscheck, and a self-inclusive value
        # would misstate replication there.
        n_cycles_total = len(cycle_tags)
        loo_denom = max(n_cycles_total - 1, 1)

        n_detected = {}
        for cyc, biochem_results in all_cycle_results:
            for r in biochem_results:
                pos = getattr(r, "position", None)
                group = getattr(r, "caap_group", "US")
                if pos is not None:
                    key = (pos, group)
                    if key not in n_detected:
                        n_detected[key] = set()
                    n_detected[key].add(cyc)

        n_detected_count = {}
        perm_pos_pval_rows = []
        for key, cycles_set in n_detected.items():
            k = len(cycles_set)
            n_detected_count[key] = k
            perm_pos_pval_rows.append({
                "Gene": gene,
                "Position": key[0],
                "caap_group": key[1],
                "n_detected": k,
                "n_cycles": n_cycles_total,
                "null_pvalue_boot": (k - 1) / loo_denom,
            })

        # ── 2. Emit raw per-(cycle, position, scheme) detail ──────────────────
        # Scoring itself is deliberately NOT done here. null_row_caas needs
        # 1 - percent_rank(null_pvalue_boot) ranked over the cycle's GENOME-WIDE
        # candidate pool, and this worker only ever sees one gene — so the rank
        # cannot be formed at this level. The parent finalizes it in pass B
        # (see _finalize_perm_scores) once every gene's rows have been counted.
        #
        # change_top/change_bottom are emitted per scheme row as 0/1 and OR-ed
        # across a position's schemes in pass B to derive change_side: a position
        # counts as "top" if ANY of its schemes changed on the top side, "bottom"
        # likewise, "both" when both hold, "none" otherwise.
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
                change_top = getattr(r, "change_top", "no_change")
                change_bottom = getattr(r, "change_bottom", "no_change")
                detail_rows.append({
                    "Gene": gene,
                    "cycle": cyc,
                    "Position": pos,
                    "caap_group": group,
                    "asr_path_score": asr_val,
                    "n_detected": n_detected_count.get((pos, group), 1),
                    "ct": 1 if change_top in _CHANGE_STATES else 0,
                    "cb": 1 if change_bottom in _CHANGE_STATES else 0,
                })

        return (gene, detail_rows, perm_pos_pval_rows)
    except Exception as e:
        logger.error(f"[perms] worker failed for {gene}: {e}", exc_info=True)
        return (gene, [], [])


def _perms_worker_wrapper(args):
    return _perms_worker(*args)


def _null_row_caas(row: Dict[str, Any], rank_lookup: Dict[str, Dict[int, float]]) -> float:
    """phen x asr for one null detail row (scoring_compute.R's caas_row).

    Single definition shared by BOTH finalize sub-passes, so the pool that
    defines the size-adjust reference and the genes scored against it cannot
    drift apart.
    """
    phen = 1.0 - rank_lookup.get(row["cycle"], {}).get(int(row["n_detected"]), 0.0)
    return phen * float(row["asr_path_score"])


def _build_cycle_score_pools(
    detail_path: Path,
    rank_lookup: Dict[str, Dict[int, float]],
) -> Dict[str, Dict[str, Any]]:
    """Sub-pass B1: per-cycle, per-direction pool of position-level null scores.

    The gene-level null statistic is size_adj_max = F(max)^n, mirroring
    scoring_compute.R section 4a. F must be the ECDF of the pool the gene's
    positions were actually drawn from, so every cycle needs its OWN pool --
    exactly as the observed side calibrates against the pool it discovered, and
    as build_percent_rank_lookup already does for the phen axis. Without this the
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

    acc: Dict[str, Dict[str, Any]] = {}

    def _bump(cyc: str, score: float, ct: int, cb: int) -> None:
        per_cycle = acc.get(cyc)
        if per_cycle is None:
            per_cycle = {k: array.array("d") for k in ("all", "top", "bottom")}
            acc[cyc] = per_cycle
        per_cycle["all"].append(score)
        if ct:
            per_cycle["top"].append(score)
        if cb:
            per_cycle["bottom"].append(score)

    def _drain(agg: Dict[Tuple[str, int], List[float]]) -> None:
        for (cyc, _pos), entry in agg.items():
            # MEAN over the position's schemes -- mirrors scoring_compute.R 2g.
            score = entry[0] / entry[1] if entry[1] else 0.0
            _bump(cyc, score, int(entry[2]), int(entry[3]))

    with gzip.open(detail_path, "rt", newline="") as f_in:
        reader = _csv.DictReader(f_in, delimiter="\t")
        current_gene: Optional[str] = None
        pos_agg: Dict[Tuple[str, int], List[float]] = {}
        for row in reader:
            gene = row["Gene"]
            if gene != current_gene:
                if current_gene is not None:
                    _drain(pos_agg)
                current_gene = gene
                pos_agg = {}
            key = (row["cycle"], int(row["Position"]))
            entry = pos_agg.setdefault(key, [0.0, 0, 0, 0])
            entry[0] += _null_row_caas(row, rank_lookup)
            entry[1] += 1
            entry[2] |= int(row["ct"])
            entry[3] |= int(row["cb"])
        if current_gene is not None:
            _drain(pos_agg)

    return {cyc: {k: array.array("d", sorted(v)) for k, v in per_cycle.items()}
            for cyc, per_cycle in acc.items()}


def _size_adj_max_null(vals: List[float], cyc: str, direction: str,
                       pools: Dict[str, Dict[str, Any]]) -> float:
    """F_cycle(max)^n -- the null mirror of scoring_compute.R's size_adj_max()."""
    if not vals:
        return 0.0
    per_cycle = pools.get(cyc)
    if per_cycle is None:
        return 0.0
    pool = per_cycle[direction]
    n_pool = len(pool)
    if n_pool == 0:
        return 0.0
    # bisect_right gives #{pool <= m}, matching findInterval() on the observed
    # side; both count ties at m as below-or-equal so the two agree exactly.
    return (bisect.bisect_right(pool, max(vals)) / n_pool) ** len(vals)


def _finalize_perm_scores(
    detail_path: Path,
    output_dir: Path,
    cycle_tags: List[str],
    rank_lookup: Dict[str, Dict[int, float]],
    sample_per_cycle_group: Optional[int] = None,
) -> None:
    """Pass B: rank within each cycle, score, and aggregate to gene x cycle stats.

    Mirrors scoring_compute.R's observed pipeline term for term:

        null_phen_score = 1 - percent_rank(null_pvalue_boot)   # within cycle i's pool
        null_row_caas   = null_phen_score * asr
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

    # Reservoir size per (cycle, scheme). Bounds both the violin sample and the
    # quantile summaries at ~K * n_cycles * 5 rows regardless of run size. The
    # full detail file stays on disk, so exact statistics remain recoverable.
    if sample_per_cycle_group is None:
        sample_per_cycle_group = int(os.environ.get("CAAS_PERMS_SAMPLE_PER_CYCLE_GROUP", "200"))

    scores_fields = ["Gene", "cycle", "global_asr", "top_asr", "bottom_asr",
                     "global_caas", "top_caas", "bottom_caas"]

    rng = random.Random(1998)
    reservoirs: Dict[Tuple[str, str], List[Tuple[str, int, float, float, float]]] = {}
    seen: Dict[Tuple[str, str], int] = {}

    # ── Sub-pass B1: per-cycle reference pools for the size-adjusted max ──────
    cycle_pools = _build_cycle_score_pools(detail_path, rank_lookup)
    logger.info("[perms] pass B1 done: size-adjust reference pools built for %d cycles",
                len(cycle_pools))

    def _q90(vals) -> float:
        return float(np.percentile(vals, 90)) if vals else 0.0

    def _flush(gene: str, pos_agg: Dict[Tuple[str, int], List[float]], writer) -> None:
        by_cycle: Dict[str, List[Tuple[float, float, str]]] = {}
        for (cyc, _pos), agg in pos_agg.items():
            asr_sum, n_schemes, caas_sum, ct, cb = agg
            asr_score = asr_sum / n_schemes if n_schemes else 0.0
            # Position score = MEAN of caas_row over the schemes that detected it
            # (scoring_compute.R section 2g); previously a 0.2-weighted sum.
            caas_score = caas_sum / n_schemes if n_schemes else 0.0
            if ct and cb:
                side = "both"
            elif ct:
                side = "top"
            elif cb:
                side = "bottom"
            else:
                side = "none"
            by_cycle.setdefault(cyc, []).append((asr_score, caas_score, side))

        rows = []
        for cyc in cycle_tags:
            items = by_cycle.get(cyc, [])
            if not items:
                # A permuted labeling that produced no hit anywhere in this gene is
                # a structural zero, not missing data -- keep the explicit 0.0 row so
                # the genes x N null matrix stays dense (scoring_caas_perms.R and the
                # report's p_zero diagnostics both rely on this).
                rows.append({"Gene": gene, "cycle": cyc,
                             "global_asr": 0.0, "top_asr": 0.0, "bottom_asr": 0.0,
                             "global_caas": 0.0, "top_caas": 0.0, "bottom_caas": 0.0})
                continue
            asr_g = [i[0] for i in items]
            asr_t = [i[0] for i in items if i[2] in ("top", "both")]
            asr_b = [i[0] for i in items if i[2] in ("bottom", "both")]
            caas_g = [i[1] for i in items]
            caas_t = [i[1] for i in items if i[2] in ("top", "both")]
            caas_b = [i[1] for i in items if i[2] in ("bottom", "both")]
            rows.append({
                "Gene": gene, "cycle": cyc,
                # ASR axis: not wired into any ranking downstream, so it stays on
                # q90 -- matching the observed gene_caas_asr_* columns, which are
                # also still q90. The two must move together if ever wired up.
                "global_asr": _q90(asr_g), "top_asr": _q90(asr_t), "bottom_asr": _q90(asr_b),
                # CAAS axis: this is what FCS actually consumes (as
                # caas_corStat_byrank), so it MUST match scoring_compute.R's
                # size_adj_max term for term or p.perm compares two different
                # statistics and silently goes wrong.
                "global_caas": _size_adj_max_null(caas_g, cyc, "all",    cycle_pools),
                "top_caas":    _size_adj_max_null(caas_t, cyc, "top",    cycle_pools),
                "bottom_caas": _size_adj_max_null(caas_b, cyc, "bottom", cycle_pools),
            })
        writer.writerows(rows)

    n_rows = 0
    with gzip.open(detail_path, "rt", newline="") as f_in, \
         open(scores_path, "w", newline="") as f_scores:
        reader = _csv.DictReader(f_in, delimiter="\t")
        writer_scores = _csv.DictWriter(f_scores, fieldnames=scores_fields, delimiter="\t")
        writer_scores.writeheader()

        current_gene: Optional[str] = None
        pos_agg: Dict[Tuple[str, int], List[float]] = {}

        for row in reader:
            gene = row["Gene"]
            if gene != current_gene:
                if current_gene is not None:
                    _flush(current_gene, pos_agg, writer_scores)
                current_gene = gene
                pos_agg = {}

            cyc = row["cycle"]
            grp = row["caap_group"]
            pos = int(row["Position"])
            asr = float(row["asr_path_score"])
            d = int(row["n_detected"])

            phen = 1.0 - rank_lookup.get(cyc, {}).get(d, 0.0)
            rc = phen * asr

            agg = pos_agg.setdefault((cyc, pos), [0.0, 0, 0.0, 0, 0])
            agg[0] += asr
            agg[1] += 1
            agg[2] += rc
            agg[3] |= int(row["ct"])
            agg[4] |= int(row["cb"])

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
            item = (gene, pos, asr, phen, rc)
            if len(res) < sample_per_cycle_group:
                res.append(item)
            else:
                j = rng.randrange(seen[key])
                if j < sample_per_cycle_group:
                    res[j] = item
            n_rows += 1

        if current_gene is not None:
            _flush(current_gene, pos_agg, writer_scores)

    # ── Sample + quantile summaries ────────────────────────────────────────────
    sample_fields = ["Gene", "Position", "caap_group", "cycle",
                     "asr_path_score", "null_phen_score", "null_row_caas"]
    with open(sample_path, "w", newline="") as f_sample:
        writer_sample = _csv.DictWriter(f_sample, fieldnames=sample_fields, delimiter="\t")
        writer_sample.writeheader()
        for (cyc, grp), res in reservoirs.items():
            writer_sample.writerows({
                "Gene": g, "Position": p, "caap_group": grp, "cycle": cyc,
                "asr_path_score": a, "null_phen_score": ph, "null_row_caas": rc,
            } for (g, p, a, ph, rc) in res)

    quant_levels = [5, 10, 25, 50, 75, 90, 95]
    quant_fields = (["cycle", "caap_group", "metric", "n_sampled", "n_total", "mean"]
                    + [f"q{q}" for q in quant_levels])
    with open(quant_path, "w", newline="") as f_quant:
        writer_quant = _csv.DictWriter(f_quant, fieldnames=quant_fields, delimiter="\t")
        writer_quant.writeheader()
        for (cyc, grp), res in reservoirs.items():
            if not res:
                continue
            for metric, idx in (("asr_path_score", 2), ("null_phen_score", 3), ("null_row_caas", 4)):
                vals = np.asarray([r[idx] for r in res], dtype=float)
                rec = {"cycle": cyc, "caap_group": grp, "metric": metric,
                       "n_sampled": len(vals), "n_total": seen.get((cyc, grp), len(vals)),
                       "mean": float(vals.mean())}
                for q, v in zip(quant_levels, np.percentile(vals, quant_levels)):
                    rec[f"q{q}"] = float(v)
                writer_quant.writerow(rec)

    logger.info(
        "[perms] pass B done: scored %d rows -> %s; sample=%s quantiles=%s",
        n_rows, scores_path.name, sample_path.name, quant_path.name,
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
    convergence_mode: str,
    workers: Optional[int],
    output_dir: Path,
    ensembl_genes_file: Optional[str] = None,
    cycles: Optional[List[str]] = None,
    max_tasks_per_child: Optional[int] = None,
) -> Path:
    """Genome-wide CAAS permulation null: load ASR once per gene, replay N permuted
    labelings, and score them the same way the observed pipeline scores itself.

    Two passes, because the null's phenotype axis is a percent_rank over each
    cycle's genome-wide candidate pool while workers only ever see a single gene:

      Pass A (parallel, one gene per worker) replays the labelings and streams raw
        per-(gene, cycle, position, scheme) detail to perm_pos_detail.tsv.gz,
        accumulating a per-cycle histogram of n_detected as rows go by.
      Pass B (single process, streaming) turns those histograms into exact
        per-cycle percent_rank lookups, re-reads the detail file, and derives
        null_row_caas -> position scores -> per-(gene, cycle) q90.

    Outputs:
      - output_dir/gene_cycle_scores.tsv   (schema unchanged; feeds caas_perms.rds)
      - output_dir/perm_pos_pval.tsv       (LOO null_pvalue_boot crosscheck table)
      - output_dir/perm_pos_detail.tsv.gz  (full per-cycle detail; re-scoring needs
                                            no ASR replay)
      - output_dir/perm_pos_quantiles.tsv  (per (cycle, scheme) distribution shape)
      - output_dir/perm_pos_sample.tsv     (cycle-stratified sample for violins)
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

    cycle_tags, cycle_trait_files, cycle_meta_files = build_cycle_inputs(
        perm_discovery_file, resample_dir, output_dir / "cycle_inputs", cycles
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

    pval_path = output_dir / "perm_pos_pval.tsv"
    detail_path = output_dir / "perm_pos_detail.tsv.gz"

    pval_fields = ["Gene", "Position", "caap_group", "n_detected", "n_cycles", "null_pvalue_boot"]
    detail_fields = ["Gene", "cycle", "Position", "caap_group", "asr_path_score", "n_detected", "ct", "cb"]

    # Generate arguments lazily
    args_generator = (
        (gene, alignment_dir, tree_file, taxid_mapping_path, asr_model,
         asr_cache_dir, posterior_threshold, convergence_mode,
         cycle_tags, cycle_trait_files, perm_discovery_file, None)
        for gene in genes
    )

    # ── Pass A: replay labelings, stream detail, count per-cycle candidate pools ─
    # hist_by_cycle[cycle][n_detected] = how many (gene, position, scheme)
    # candidates that cycle discovered at that replication level. Bounded by
    # n_cycles^2 counters (~1000x1000 ints, a few MB) no matter how many rows the
    # detail file holds.
    hist_by_cycle: Dict[str, Dict[int, int]] = {}
    pool = mp.Pool(processes=effective_workers, maxtasksperchild=maxtasks)
    n_genes = 0
    n_detail_rows = 0
    try:
        results_iterator = pool.imap_unordered(_perms_worker_wrapper, args_generator, chunksize=10)

        with open(pval_path, "w", newline="") as f_pval, \
             gzip.open(detail_path, "wt", newline="") as f_detail:

            writer_pval = _csv.DictWriter(f_pval, fieldnames=pval_fields, delimiter="\t")
            writer_detail = _csv.DictWriter(f_detail, fieldnames=detail_fields, delimiter="\t")

            writer_pval.writeheader()
            writer_detail.writeheader()

            for _gene, detail_rows, pval_rows in results_iterator:
                if not detail_rows:
                    continue
                writer_pval.writerows(pval_rows)
                writer_detail.writerows(detail_rows)
                n_detail_rows += len(detail_rows)
                for row in detail_rows:
                    cyc_hist = hist_by_cycle.setdefault(row["cycle"], {})
                    d = row["n_detected"]
                    cyc_hist[d] = cyc_hist.get(d, 0) + 1
                n_genes += 1
    finally:
        pool.close()
        pool.join()

    logger.info(
        f"[perms] pass A done: {n_genes} genes, {n_detail_rows} (gene,cycle,position,scheme) "
        f"rows across {len(hist_by_cycle)} cycles -> {detail_path.name}"
    )
    if n_genes < len(genes):
        logger.warning(
            f"[perms] pass A scored {n_genes}/{len(genes)} genes with per-cycle CAAS — "
            f"{len(genes) - n_genes} contributed nothing (no cached ASR, or no CAAS "
            f"survived any replayed labeling). If asr_mode=compute, ensure the live "
            f"disambiguation completed before this step and covers every gene that "
            f"appears in the permulation discovery."
        )
    if n_detail_rows == 0:
        logger.error(
            "[perms] pass A produced ZERO detail rows — the permulation null is empty. "
            "Downstream caas_perms.rds / FCS p.perm will be degenerate. Most common "
            "cause: asr_mode=compute with the ASR cache not yet populated when this "
            "step ran (see the gate in caas_permulation.nf / CAAS_PERMULATION)."
        )

    # ── Pass B: rank within each cycle, score, aggregate ────────────────────────
    rank_lookup = build_percent_rank_lookup(hist_by_cycle)
    pool_sizes = {c: sum(h.values()) for c, h in hist_by_cycle.items()}
    if pool_sizes:
        sizes = sorted(pool_sizes.values())
        logger.info(
            "[perms] per-cycle candidate pool size: min=%d median=%d max=%d "
            "(observed side ranks over its own discovered pool the same way)",
            sizes[0], sizes[len(sizes) // 2], sizes[-1],
        )

    _finalize_perm_scores(
        detail_path=detail_path,
        output_dir=output_dir,
        cycle_tags=cycle_tags,
        rank_lookup=rank_lookup,
    )

    logger.info(f"[perms] successfully aggregated {n_genes} genes to summaries inside {output_dir}")
    return output_dir


# Backward-compatible alias
convert_biochem_result_to_dict = convert_convergence_result_to_dict
