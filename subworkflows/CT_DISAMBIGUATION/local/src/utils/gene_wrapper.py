#!/usr/bin/env python3
"""
Gene Processing Wrapper with Multiprocessing

Orchestrates parallel processing of genes using existing single_gene_pipeline
logic. Streams lightweight per-position result dicts into an aggregation SQLite DB,
then exports CSV/JSON from the DB.

Author: ASR Integration
Date: 2025-12-03 (revised 2025-12-09)
"""

import logging
import multiprocessing as mp
import os
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
            return (gene, [], [], [])

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
            return (gene, [], [], [])

        # ── 1. Calculate recovery p-values on the fly ──────────────────────────
        scheme_weights_int = {"US": 0.2, "GS4": 0.2, "GS3": 0.2, "GS2": 0.2, "GS1": 0.2}
        n_cycles_total = len(cycle_tags)
        
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

        recovery_pval = {}
        perm_pos_pval_rows = []
        for key, cycles_set in n_detected.items():
            k = len(cycles_set)
            pval = (k + 1) / (n_cycles_total + 1)
            recovery_pval[key] = pval
            perm_pos_pval_rows.append({
                "Gene": gene,
                "Position": key[0],
                "caap_group": key[1],
                "n_detected": k,
                "n_cycles": n_cycles_total,
                "recovery_pval": pval
            })

        # ── 2. Group by (cycle, Position) and calculate asr+caas scores ────────
        pos_groups = {}
        for cyc, biochem_results in all_cycle_results:
            for r in biochem_results:
                pos = getattr(r, "position", None)
                if pos is not None:
                    key = (cyc, pos)
                    if key not in pos_groups:
                        pos_groups[key] = []
                    pos_groups[key].append(r)

        pos_scores_by_cycle = {}
        for (cyc, pos), r_list in pos_groups.items():
            asr_vals = []
            row_caas_sum = 0.0
            has_change_top = False
            has_change_bottom = False
            
            for r in r_list:
                group = getattr(r, "caap_group", "US")
                asr_val = getattr(r, "asr_path_score", 0.0)
                if asr_val is None:
                    asr_val = 0.0
                change_top = getattr(r, "change_top", "no_change")
                change_bottom = getattr(r, "change_bottom", "no_change")
                
                pval = recovery_pval[(pos, group)]
                weight = scheme_weights_int.get(group, 0.2)
                
                asr_vals.append(asr_val)
                row_caas_sum += (1.0 - pval) * asr_val * weight
                
                if change_top in ("convergent", "codivergent", "divergent"):
                    has_change_top = True
                if change_bottom in ("convergent", "codivergent", "divergent"):
                    has_change_bottom = True

            asr_score = sum(asr_vals) / len(asr_vals) if asr_vals else 0.0
            
            if has_change_top and has_change_bottom:
                change_side = "both"
            elif has_change_top:
                change_side = "top"
            elif has_change_bottom:
                change_side = "bottom"
            else:
                change_side = "none"

            if cyc not in pos_scores_by_cycle:
                pos_scores_by_cycle[cyc] = []
            pos_scores_by_cycle[cyc].append((asr_score, row_caas_sum, change_side))

        # ── 3. Calculate 90th percentile scores per cycle (gene_cycle_rows) ────
        import numpy as np
        gene_cycle_rows = []
        for cyc in cycle_tags:
            items = pos_scores_by_cycle.get(cyc, [])
            if not items:
                gene_cycle_rows.append({
                    "Gene": gene,
                    "cycle": cyc,
                    "global_asr": 0.0, "top_asr": 0.0, "bottom_asr": 0.0,
                    "global_caas": 0.0, "top_caas": 0.0, "bottom_caas": 0.0
                })
                continue
                
            asr_vals_global = [item[0] for item in items]
            asr_vals_top = [item[0] for item in items if item[2] in ("top", "both")]
            asr_vals_bottom = [item[0] for item in items if item[2] in ("bottom", "both")]
            
            caas_vals_global = [item[1] for item in items]
            caas_vals_top = [item[1] for item in items if item[2] in ("top", "both")]
            caas_vals_bottom = [item[1] for item in items if item[2] in ("bottom", "both")]
            
            def q90(vals):
                return float(np.percentile(vals, 90)) if vals else 0.0

            gene_cycle_rows.append({
                "Gene": gene,
                "cycle": cyc,
                "global_asr": q90(asr_vals_global),
                "top_asr": q90(asr_vals_top),
                "bottom_asr": q90(asr_vals_bottom),
                "global_caas": q90(caas_vals_global),
                "top_caas": q90(caas_vals_top),
                "bottom_caas": q90(caas_vals_bottom)
            })

        # ── 4. Sub-sample rows for distribution plotting ───────────────────────
        import random
        # Seed locally inside the worker for deterministic sub-sampling
        random.seed(hash(gene) & 0xffffffff)
        
        perm_pos_sample_rows = []
        by_group = {g: [] for g in ["US", "GS1", "GS2", "GS3", "GS4"]}
        for cyc, biochem_results in all_cycle_results:
            for r in biochem_results:
                pos = getattr(r, "position", None)
                if pos is not None:
                    group = getattr(r, "caap_group", "US")
                    asr_val = getattr(r, "asr_path_score", 0.0)
                    if asr_val is None:
                        asr_val = 0.0
                    pval = recovery_pval[(pos, group)]
                    weight = scheme_weights_int.get(group, 0.2)
                    row_caas = (1.0 - pval) * asr_val * weight
                    
                    by_group[group].append({
                        "Gene": gene,
                        "Position": pos,
                        "caap_group": group,
                        "cycle": cyc,
                        "asr_path_score": asr_val,
                        "recovery_pval": pval,
                        "row_caas": row_caas
                    })
        for g, rows in by_group.items():
            if len(rows) > 5:
                perm_pos_sample_rows.extend(random.sample(rows, 5))
            else:
                perm_pos_sample_rows.extend(rows)

        return (gene, gene_cycle_rows, perm_pos_pval_rows, perm_pos_sample_rows)
    except Exception as e:
        logger.error(f"[perms] worker failed for {gene}: {e}", exc_info=True)
        return (gene, [], [], [])


def _perms_worker_wrapper(args):
    return _perms_worker(*args)


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
    labelings, and stream on-the-fly aggregated summaries to output TSV files.

    Outputs:
      - output_dir/gene_cycle_scores.tsv
      - output_dir/perm_pos_pval.tsv
      - output_dir/perm_pos_sample.tsv
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

    scores_path = output_dir / "gene_cycle_scores.tsv"
    pval_path = output_dir / "perm_pos_pval.tsv"
    sample_path = output_dir / "perm_pos_sample.tsv"

    scores_fields = ["Gene", "cycle", "global_asr", "top_asr", "bottom_asr", "global_caas", "top_caas", "bottom_caas"]
    pval_fields = ["Gene", "Position", "caap_group", "n_detected", "n_cycles", "recovery_pval"]
    sample_fields = ["Gene", "Position", "caap_group", "cycle", "asr_path_score", "recovery_pval", "row_caas"]

    # Generate arguments lazily
    args_generator = (
        (gene, alignment_dir, tree_file, taxid_mapping_path, asr_model,
         asr_cache_dir, posterior_threshold, convergence_mode,
         cycle_tags, cycle_trait_files, perm_discovery_file, None)
        for gene in genes
    )

    pool = mp.Pool(processes=effective_workers, maxtasksperchild=maxtasks)
    n_genes = 0
    try:
        results_iterator = pool.imap_unordered(_perms_worker_wrapper, args_generator, chunksize=10)
        
        with open(scores_path, "w", newline="") as f_scores, \
             open(pval_path, "w", newline="") as f_pval, \
             open(sample_path, "w", newline="") as f_sample:
             
             writer_scores = _csv.DictWriter(f_scores, fieldnames=scores_fields, delimiter="\t")
             writer_pval = _csv.DictWriter(f_pval, fieldnames=pval_fields, delimiter="\t")
             writer_sample = _csv.DictWriter(f_sample, fieldnames=sample_fields, delimiter="\t")
             
             writer_scores.writeheader()
             writer_pval.writeheader()
             writer_sample.writeheader()
             
             for _gene, scores_rows, pval_rows, sample_rows in results_iterator:
                 if scores_rows:
                     writer_scores.writerows(scores_rows)
                     writer_pval.writerows(pval_rows)
                     writer_sample.writerows(sample_rows)
                     n_genes += 1
    finally:
        pool.close()
        pool.join()

    logger.info(f"[perms] successfully aggregated {n_genes} genes to summaries inside {output_dir}")
    return output_dir


# Backward-compatible alias
convert_biochem_result_to_dict = convert_convergence_result_to_dict
