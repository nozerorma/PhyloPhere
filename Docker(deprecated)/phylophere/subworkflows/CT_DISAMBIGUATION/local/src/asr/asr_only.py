#!/usr/bin/env python3
"""
ASR-only runner for single-gene analyses.

Runs ancestral state reconstruction for one or more genes without triggering
downstream biochemical/reporting steps. Useful to precompute ASR results that
the main pipeline can later reuse with --asr-mode precomputed.
"""

from pathlib import Path
import argparse
import logging
import json
import sys

from src.asr.asr_single import (
    run_asr_pipeline,
    SingleGeneASRConfig,
    validate_asr_inputs,
)
from src.utils.io_utils import find_gene_alignment

# Add src to path for utility imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))
from src.utils.logger import configure_logging

logger = logging.getLogger(__name__)


def _run_gene_asr(
    gene: str,
    alignment_dir: Path,
    tree: Path,
    taxid_mapping: Path,
    output_dir: Path,
    asr_model: str,
    threads: int,
    posterior_threshold: float,
    run_diagnostics: bool,
    skip_if_exists: bool,
    ensembl_genes,
):
    """Worker-friendly wrapper to run ASR for a single gene."""
    alignment_path = find_gene_alignment(alignment_dir, gene, ensembl_genes)

    gene_out = output_dir
    gene_out.mkdir(parents=True, exist_ok=True)

    config = SingleGeneASRConfig(
        alignment_path=alignment_path,
        tree_path=tree,
        taxid_path=taxid_mapping,
        output_dir=gene_out,
        model=asr_model,
        threads=threads,
        posterior_threshold=posterior_threshold,
        run_diagnostics=run_diagnostics,
    )

    validate_asr_inputs(config)

    logger.info(f"Running ASR for {gene} -> {gene_out}")
    from src.utils.concurrency import codeml_slot

    with codeml_slot():
        result = run_asr_pipeline(
            gene=gene, config=config, skip_if_exists=skip_if_exists
        )
    return {
        "gene": gene,
        "success": bool(result),
        "rst_file": str(result.rst_file) if result and result.rst_file else None,
        "tree_file": str(result.tree_file) if result and result.tree_file else None,
        "node_id_map": (
            str(result.node_id_map) if result and result.node_id_map else None
        ),
        "output_dir": str(gene_out),
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Run ASR only for one or more genes")
    parser.add_argument(
        "--genes", "-g", nargs="+", required=True, help="Gene names to process"
    )
    parser.add_argument(
        "--alignment_dir",
        "-a",
        type=Path,
        required=True,
        help="Directory with alignments",
    )
    parser.add_argument(
        "--ensembl-genes-file",
        type=Path,
        help="TSV/CSV with gene column to enforce exact alignment matching",
    )
    parser.add_argument(
        "--tree", "-t", type=Path, required=True, help="Tree file (.nex/.nwk)"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        default=Path("./asr_only_output"),
        help="Output directory",
    )
    parser.add_argument("--taxid_mapping", type=Path, help="TaxID mapping file")
    parser.add_argument("--asr_model", default="lg", help="ASR model (default: lg)")
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads for ASR (default: 1)"
    )
    parser.add_argument(
        "--workers",
        "-w",
        type=int,
        default=None,
        help="Parallel ASR workers across genes (default: auto)",
    )
    parser.add_argument(
        "--codeml-concurrency",
        type=int,
        default=None,
        help="Max concurrent codeml runs (default: auto)",
    )
    parser.add_argument(
        "--posterior-threshold",
        type=float,
        default=0.7,
        help="Posterior threshold (default: 0.7)",
    )
    parser.add_argument(
        "--skip-if-exists",
        action="store_true",
        help="Skip ASR if outputs already present",
    )
    parser.add_argument(
        "--run-diagnostics",
        action="store_true",
        help="Enable diagnostics: node-level posteriors and per-position tip details in diagnostics/",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional path to write ASR logs (overwrites if exists)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    configure_logging(verbose=args.verbose, log_file=args.log_file)

    if not args.alignment_dir.exists():
        logger.error(f"Alignment dir not found: {args.alignment_dir}")
        sys.exit(1)
    if not args.tree.exists():
        logger.error(f"Tree file not found: {args.tree}")
        sys.exit(1)
    if args.taxid_mapping and not args.taxid_mapping.exists():
        logger.error(f"TaxID mapping file not found: {args.taxid_mapping}")
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []

    ensembl_genes = None
    if args.ensembl_genes_file:
        try:
            from src.data.loaders import load_ensembl_genes

            ensembl_genes = load_ensembl_genes(args.ensembl_genes_file)
            logger.info(
                f"Loaded {len(ensembl_genes)} genes from {args.ensembl_genes_file}"
            )
        except Exception as exc:
            logger.error(f"Failed to load Ensembl genes file: {exc}")
            sys.exit(1)

    from concurrent.futures import ProcessPoolExecutor, as_completed
    from src.utils.concurrency import plan_concurrency, init_worker
    import multiprocessing as mp

    effective_workers, threads = plan_concurrency(args.workers, args.threads, logger)

    if effective_workers > 1:
        manager = mp.Manager()
        codeml_sem = (
            manager.Semaphore(max(1, args.codeml_concurrency))
            if args.codeml_concurrency
            else None
        )
        with ProcessPoolExecutor(
            max_workers=effective_workers,
            initializer=init_worker,
            initargs=(threads, codeml_sem),
        ) as executor:
            future_to_gene = {
                executor.submit(
                    _run_gene_asr,
                    gene,
                    args.alignment_dir,
                    args.tree,
                    args.taxid_mapping,
                    args.output_dir,
                    args.asr_model,
                    threads,
                    args.posterior_threshold,
                    args.run_diagnostics,
                    args.skip_if_exists,
                    ensembl_genes,
                ): gene
                for gene in args.genes
            }
            for fut in as_completed(future_to_gene):
                gene = future_to_gene[fut]
                try:
                    summaries.append(fut.result())
                except FileNotFoundError as e:
                    logger.error(f"{e}")
                except Exception as exc:
                    logger.error(f"ASR failed for {gene}: {exc}", exc_info=True)
    else:
        for gene in args.genes:
            try:
                summaries.append(
                    _run_gene_asr(
                        gene,
                        args.alignment_dir,
                        args.tree,
                        args.taxid_mapping,
                        args.output_dir,
                        args.asr_model,
                        threads,
                        args.posterior_threshold,
                        args.run_diagnostics,
                        args.skip_if_exists,
                        ensembl_genes,
                    )
                )
            except FileNotFoundError as e:
                logger.error(f"{e}")
            except Exception as exc:
                logger.error(f"ASR failed for {gene}: {exc}", exc_info=True)

    summary_path = args.output_dir / "asr_only_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2))
    logger.info(f"Wrote ASR-only summary to {summary_path}")


if __name__ == "__main__":
    main()
