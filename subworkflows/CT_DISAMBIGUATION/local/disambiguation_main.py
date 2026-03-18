#!/usr/bin/env python3
"""CLI entry point for the CAAS aggregation/disambiguation pipeline."""

import sys
import argparse
import logging
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.data.loaders import read_caas_metadata_table
from src.utils.gene_wrapper import process_all_genes
from src.reporting.disambiguation_writers import (
    write_caas_convergence_csvs,
)
from src.reporting.disambiguation_json import (
    export_aggregated_convergence_json,
    export_gene_summaries_json,
)
from src.reporting.gene_lists import export_gene_lists
from src.plots.bulk_plots import generate_bulk_plots
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse CLI arguments for the aggregation pipeline.

    :returns: Parsed command-line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description="CAAS Aggregation Pipeline - Process all genes with ASR-based convergence analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required inputs
    parser.add_argument(
        "--alignment-dir", required=True, help="Directory containing alignment files"
    )
    parser.add_argument(
        "--tree", required=True, help="Phylogenetic tree file (Newick/Nexus)"
    )
    parser.add_argument(
        "--caas-metadata", required=True, help="CAAS metadata file (all genes)"
    )
    parser.add_argument(
        "--trait-file",
        required=True,
        help="Trait file with species/contrast/trait/pair columns",
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for results"
    )
    parser.add_argument(
        "--ensembl-genes-file",
        default=None,
        help="TSV/CSV with gene column; used to strictly match alignment prefixes",
    )

    # ASR configuration
    parser.add_argument(
        "--asr-mode",
        choices=["compute", "precomputed"],
        default="precomputed",
        help="ASR mode (default: precomputed)",
    )
    parser.add_argument(
        "--asr-model", default="lg", help="Substitution model for ASR (default: lg)"
    )
    parser.add_argument(
        "--asr-cache-dir", default=None, help="Directory for precomputed ASR files"
    )
    parser.add_argument(
        "--posterior-threshold",
        type=float,
        default=0.0,
        help="Posterior probability threshold (default: 0.0)",
    )

    # Analysis configuration
    parser.add_argument(
        "--convergence-mode",
        choices=["focal_clade", "mrca"],
        default="mrca",
        help="Convergence analysis mode (default: mrca)",
    )
    parser.add_argument(
        "--taxid-mapping", default=None, help="TaxID to species mapping file"
    )

    # Performance
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: cpu_count-2)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for codeml per gene (default: 1)",
    )
    parser.add_argument(
        "--codeml-concurrency",
        type=int,
        default=None,
        help="Max concurrent codeml/ASR runs (default: auto based on workers and threads)",
    )
    parser.add_argument(
        "--max-tasks-per-child",
        type=int,
        default=None,
        help="Max tasks per worker process before restart (maxtasksperchild). Overrides CAAS_MAX_TASKS_PER_CHILD env var if set.",
    )

    # Filters and options
    parser.add_argument(
        "--include-non-significant",
        action="store_true",
        help="Include non-significant CAAS in outputs",
    )
    parser.add_argument(
        "--skip-gene-lists",
        action="store_true",
        help="Skip per-pattern gene list export",
    )
    parser.add_argument(
        "--allow-low-confidence",
        action="store_true",
        help="Continue with warnings on low ASR confidence",
    )
    # Diagnostics output
    parser.add_argument(
        "--run-diagnostics",
        action="store_true",
        help="Enable diagnostics output: node-level posteriors, tip residue details to diagnostics/ subdirectory",
    )

    # Logging
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional path to write aggregated logs (overwrites if exists)",
    )

    return parser.parse_args()


def _compute_max_pairs_from_trait(trait_file: Path) -> int:
    import csv

    max_pair_id = 0
    with open(trait_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pair_val = row.get("pair")
            # Skip missing or empty values before attempting conversion
            if pair_val is None or str(pair_val).strip() == "":
                continue
            try:
                pair_id = int(str(pair_val).strip())
                max_pair_id = max(max_pair_id, pair_id)
            except Exception:
                continue
    return max_pair_id or 1


def main():
    """Run the aggregation pipeline end-to-end."""
    args = parse_arguments()

    configure_logging(verbose=args.verbose, log_file=args.log_file)

    logger.info("=" * 80)
    logger.info("CAAS Aggregation Pipeline")
    logger.info("=" * 80)

    # Log configuration
    logger.info("Configuration:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Reduce noise from verbose submodules
    for noisy in [
        # "src.utils.gene_wrapper",
        "src.phylo.species_mapping",
    ]:
        logging.getLogger(noisy).setLevel(logging.WARNING)

    # Load metadata to discover genes
    logger.info(f"Loading CAAS metadata from {args.caas_metadata}")
    start_time = time.time()

    try:
        metadata_df = read_caas_metadata_table(Path(args.caas_metadata))
        unique_genes = (
            metadata_df["GenePos"]
            .apply(lambda x: x.split("_")[0] if "_" in str(x) else None)
            .dropna()
            .unique()
            .tolist()
        )

        logger.info(f"Found {len(unique_genes)} unique genes in metadata")
        logger.info(f"Total CAAS positions: {len(metadata_df)}")

        if not args.include_non_significant:
            significant_count = metadata_df["isSignificant"].sum()
            logger.info(f"Significant CAAS: {significant_count}")

    except Exception as e:
        logger.error(f"Failed to load metadata: {e}")
        raise

    load_time = time.time() - start_time
    logger.info(f"Metadata loaded in {load_time:.2f}s")

    # Determine schema max_pairs from trait file (fixed schema requirement)
    try:
        max_pairs = _compute_max_pairs_from_trait(Path(args.trait_file))
        logger.info(f"Detected max_pairs={max_pairs} from trait file")
    except Exception as e:
        logger.error(f"Failed to compute max_pairs from trait file: {e}")
        max_pairs = None

    # Process all genes
    logger.info(f"Processing {len(unique_genes)} genes...")
    process_start = time.time()

    try:
        proc_res = process_all_genes(
            genes=unique_genes,
            alignment_dir=args.alignment_dir,
            tree_file=args.tree,
            caas_metadata_path=args.caas_metadata,
            trait_file_path=args.trait_file,
            taxid_mapping_path=args.taxid_mapping,
            asr_mode=args.asr_mode,
            asr_model=args.asr_model,
            asr_cache_dir=args.asr_cache_dir,
            posterior_threshold=args.posterior_threshold,
            convergence_mode=args.convergence_mode,
            threads_per_gene=args.threads,
            workers=args.workers,
            max_tasks_per_child=args.max_tasks_per_child,
            include_non_significant=args.include_non_significant,
            run_diagnostics=args.run_diagnostics,
            output_dir=output_dir,
            ensembl_genes_file=args.ensembl_genes_file,
            max_codeml=args.codeml_concurrency,
        )
        # process_all_genes now returns (caas_results, export_info)
        if isinstance(proc_res, tuple) and len(proc_res) == 2:
            caas_results, export_info = proc_res
        else:
            caas_results = proc_res
            export_info = None
    except Exception as e:
        logger.error(f"Gene processing failed: {e}")
        raise

    process_time = time.time() - process_start
    logger.info(f"Gene processing completed in {process_time:.2f}s")
    logger.info(f"  CAAS results: {len(caas_results)}")
    logger.info("  Conserved-pair ASR flags carried in master CSV rows")

    # Write output CSVs
    logger.info("Writing output files...")
    write_start = time.time()
    try:
        # If process returned export_info (DB-backed), skip writing CSVs again
        if export_info:
            caas_files = [Path(p) for p in export_info.get("caas_files", [])]
            logger.info(
                f"  Skipping CSV generation; files produced by DB exporter: {caas_files}"
            )
            # If exporter produced a summary JSON, log it
            summary_json = export_info.get("summary_json")
            if summary_json:
                logger.info(f"  Aggregated JSON summary: {summary_json}")
        else:
            caas_files = write_caas_convergence_csvs(
                results=caas_results,
                output_dir=output_dir,
                max_pairs=max_pairs,
            )
            logger.info(f"  CAAS convergence CSVs: {len(caas_files)} files")
            # Export aggregated convergence JSON
            try:
                if caas_results:
                    json_file = export_aggregated_convergence_json(
                        caas_results=caas_results,
                        output_path=output_dir / "caas_convergence_summary.json",
                    )
                    logger.info(f"  Convergence JSON: {json_file}")

                    # Per-gene JSON summaries
                    json_dir = output_dir / "json_summaries"
                    gene_jsons = export_gene_summaries_json(
                        caas_results=caas_results,
                        output_dir=json_dir,
                    )
                    logger.info(
                        f"  Per-gene JSONs: {len(gene_jsons)} files in {json_dir}"
                    )
            except Exception as e:
                logger.warning(f"JSON export skipped: {e}")

        # Export gene lists per pattern (all and significant subsets)
        if not args.skip_gene_lists:
            try:
                # If process returned in-memory CAAS results, use them. Otherwise
                # attempt to load the master CSV produced by the DB exporter so
                # gene lists are generated even in DB-backed mode.
                if caas_results:
                    source_results = caas_results
                else:
                    source_results = []
                    try:
                        import csv as _csv

                        master_csv = (
                            caas_files[0]
                            if caas_files
                            else (output_dir / "caas_convergence_master.csv")
                        )
                        if Path(master_csv).exists():
                            with open(master_csv, "r", encoding="utf-8") as _f:
                                reader = _csv.DictReader(_f)
                                for row in reader:
                                    # Normalize common boolean-like fields
                                    for bool_field in (
                                        "is_significant",
                                        "is_stable",
                                        "isAntiCAAS",
                                        "ASR_ConsAntiCAAS",
                                    ):
                                        if bool_field in row:
                                            val = row[bool_field]
                                            if isinstance(val, str):
                                                v = val.strip().lower()
                                                row[bool_field] = v in (
                                                    "true",
                                                    "1",
                                                    "yes",
                                                    "y",
                                                )
                                            else:
                                                row[bool_field] = bool(val)
                                    # Normalize position fields: CSV uses `msa_pos`; gene lists expect `position_zero_based` or `position`.
                                    msa_pos = row.get("msa_pos") or row.get("position")
                                    if msa_pos is not None and msa_pos != "":
                                        # keep as string/int as-is; downstream code stringifies
                                        row["position_zero_based"] = msa_pos
                                    else:
                                        # Ensure key exists to avoid missing field issues
                                        row.setdefault("position_zero_based", "")

                                    # Ensure pattern_type exists (some CSV exports may call it 'pattern')
                                    if "pattern_type" not in row and "pattern" in row:
                                        row["pattern_type"] = row.get("pattern")

                                    source_results.append(row)
                    except Exception as e:
                        logger.warning(
                            f"Failed to read master CSV for gene list export: {e}"
                        )

                if source_results:
                    gene_list_paths = export_gene_lists(
                        results=source_results, output_root=output_dir
                    )
                    all_pattern_count = len(
                        gene_list_paths.get("all", {}).get("by_pattern", {})
                    )
                    sig_pattern_count = len(
                        gene_list_paths.get("significant", {}).get("by_pattern", {})
                    )
                    all_change_count = len(
                        gene_list_paths.get("all", {}).get("by_change", {})
                    )
                    sig_change_count = len(
                        gene_list_paths.get("significant", {}).get("by_change", {})
                    )
                    logger.info(
                        f"  Gene lists: by_pattern (all={all_pattern_count}, sig={sig_pattern_count}), "
                        f"by_change (all={all_change_count}, sig={sig_change_count})"
                    )
            except Exception as e:
                logger.warning(f"Gene list export skipped: {e}")

        asr_root = (
            Path(args.asr_cache_dir) if args.asr_cache_dir else output_dir / "asr"
        )
        node_dumps_root = output_dir / "diagnostics" / "node_dumps"

        # Bulk plots (best-effort; skip if matplotlib unavailable). Generate even if
        # CAAS list is empty as long as the master CSV exists, so the plots directory
        # is created and we emit a clear log message when there is no data.
        try:
            master_csv_raw = (
                caas_files[0]
                if caas_files
                else (output_dir / "caas_convergence_master.csv")
            )
            master_csv = Path(master_csv_raw)
            ensembl_path = (
                Path(args.ensembl_genes_file) if args.ensembl_genes_file else None
            )

            if master_csv.exists():
                generate_bulk_plots(
                    caas_csv=master_csv,
                    output_dir=output_dir / "plots_bulk",
                    ensembl_csv=ensembl_path,
                    include_non_significant=args.include_non_significant,
                    asr_root=asr_root,
                    node_dumps_root=node_dumps_root,
                )
            else:
                logger.info(f"Skipping bulk plots; master CSV not found: {master_csv}")
        except Exception as e:
            logger.warning(f"Bulk plot generation skipped: {e}")

    except Exception as e:
        logger.error(f"Failed to write outputs: {e}")
        raise

    write_time = time.time() - write_start
    logger.info(f"Outputs written in {write_time:.2f}s")

    # Summary
    total_time = time.time() - start_time
    logger.info("=" * 80)
    logger.info("Pipeline completed successfully!")
    logger.info(f"Total time: {total_time:.2f}s")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
