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
from src.data.loaders import list_gene_caas_positions
from src.utils.io_utils import find_gene_alignment

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
            if not hasattr(ns, "position_zero_based"):
                if "msa_pos" in result and result["msa_pos"] is not None:
                    ns.position_zero_based = result["msa_pos"]
                elif "position" in result and result["position"] is not None:
                    try:
                        ns.position_zero_based = int(result["position"]) - 1
                    except Exception:
                        ns.position_zero_based = result["position"]

            if not hasattr(ns, "position") and "position" in result:
                ns.position = result.get("position")

            if not hasattr(ns, "caas_pvalue") and "pvalue" in result:
                ns.caas_pvalue = result.get("pvalue")

            if not hasattr(ns, "pvalue_boot") and "pvalue_boot" in result:
                ns.pvalue_boot = result.get("pvalue_boot")
            if not hasattr(ns, "pvalue_boot") and "pvalue.boot" in result:
                ns.pvalue_boot = result.get("pvalue.boot")

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

            # Low confidence nodes extraction
            if not hasattr(ns, "low_confidence_nodes"):
                nsd = result.get("node_state_details") or {}
                if isinstance(nsd, dict) and "low_confidence_nodes" in nsd:
                    ns.low_confidence_nodes = nsd.get("low_confidence_nodes")
                else:
                    ns.low_confidence_nodes = result.get("low_confidence_nodes")

            result = ns
        except Exception:
            # If normalization fails, keep original; downstream getattr will be defensive
            pass

    # Core identity
    result_dict: Dict[str, Any] = {
        "gene": getattr(result, "gene", None),
        "msa_pos": getattr(result, "position_zero_based", None),  # 0-based
        "position": getattr(result, "position", None),  # usually 1-based
        "tag": getattr(result, "tag", None),
        "caas": getattr(result, "caas", None),
        "is_significant": bool(getattr(result, "is_significant", False)),
        "pvalue": getattr(result, "caas_pvalue", None),
        "pvalue_boot": getattr(result, "pvalue_boot", None),
        "caap_group": getattr(result, "caap_group", "US"),
        "amino_encoded": getattr(result, "amino_encoded", ""),
        "is_conserved_meta": bool(getattr(result, "is_conserved_meta", False)),
        "conserved_pair": getattr(result, "conserved_pair", ""),
        "sig_hyp": getattr(result, "sig_hyp", None),
        "sig_perm": getattr(result, "sig_perm", None),
        "sig_both": getattr(result, "sig_both", None),
        "multi_hypothesis": multi_hypothesis,
        "comments": "",
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
    result_dict["pattern_type"] = getattr(result, "pattern_type", None)

    # Change tracking
    result_dict["change_top"] = getattr(result, "change_top", "no_change")
    result_dict["change_bottom"] = getattr(result, "change_bottom", "no_change")
    result_dict["change_side"] = getattr(result, "change_side", "none")
    result_dict["parallel_top"] = getattr(result, "parallel_top", None)
    result_dict["parallel_bottom"] = getattr(result, "parallel_bottom", None)
    result_dict["parallel_type"] = getattr(result, "parallel_type", "none")

    # Low confidence nodes
    lcn = getattr(result, "low_confidence_nodes", None)
    if lcn:
        if isinstance(lcn, (list, tuple, set)):
            result_dict["low_confidence_nodes"] = ",".join(str(x) for x in lcn)
        else:
            result_dict["low_confidence_nodes"] = str(lcn)

    # Conserved-pair validation flags (set upstream)
    result_dict["asr_is_conserved"] = bool(getattr(result, "asr_is_conserved", False))
    result_dict["asr_root_conserved"] = bool(
        getattr(result, "asr_root_conserved", False)
    )

    # Keep rich payloads for downstream viz/reporting
    if getattr(result, "pair_details", None):
        result_dict["pair_details"] = getattr(result, "pair_details")
    if getattr(result, "node_mapping", None):
        result_dict["node_mapping"] = getattr(result, "node_mapping")
    if getattr(result, "node_state_details", None):
        result_dict["node_state_details"] = getattr(result, "node_state_details")

    result_dict["ambiguous"] = bool(getattr(result, "ambiguous", False))

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

    # Merge low confidence nodes
    all_lcn = []
    for r in results_group:
        low_conf_nodes = getattr(r, "low_confidence_nodes", None)
        if low_conf_nodes:
            if isinstance(low_conf_nodes, (list, tuple, set)):
                all_lcn.extend(list(low_conf_nodes))
            else:
                all_lcn.append(low_conf_nodes)
    if all_lcn:
        merged["low_confidence_nodes"] = ",".join(sorted(set(str(x) for x in all_lcn)))

    # Merge CAAS labels
    caass = [
        getattr(r, "caas", None) for r in results_group if getattr(r, "caas", None)
    ]
    if caass:
        merged["caas_merged"] = "; ".join(sorted(set(str(x) for x in caass)))

    # Merge boolean/meta flags across hypothesis rows
    merged["asr_is_conserved"] = any(
        bool(getattr(r, "asr_is_conserved", False)) for r in results_group
    )
    merged["asr_root_conserved"] = any(
        bool(getattr(r, "asr_root_conserved", False)) for r in results_group
    )
    merged["is_conserved_meta"] = any(
        bool(getattr(r, "is_conserved_meta", False)) for r in results_group
    )

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
    include_non_significant: bool,
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

        logger.info(f"Processing gene: {gene} ({len(caas_positions)} CAAS positions)")

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
                logger.info(
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
                    logger.info(f"Wrote posterior JSONL to {posterior_dump_jsonl}")
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
            include_non_significant=include_non_significant,
        )

        if not include_non_significant:
            biochem_results = [
                r for r in biochem_results if getattr(r, "is_significant", False)
            ]

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
                    msa_pos = getattr(r, "position_zero_based", None)
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

        logger.info(f"  {gene}: {len(biochem_results)} CAAS processed (streamed to DB)")

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
    include_non_significant: bool,
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
                        include_non_significant,
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
        skipped_file = output_dir / "skipped_genes.txt"
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


# Backward-compatible alias
convert_biochem_result_to_dict = convert_convergence_result_to_dict
