#!/usr/bin/env python3
"""
CT_ACCUMULATION entry point — adapted for PhyloPhere CT_ACCUMULATION module.

Orchestrates two phases:
  1. aggregate  — build global position matrix and per-group position CSVs
  2. randomize  — permutation test for per-gene CAAS accumulation

Changes vs. original:
  - Randomize phase now exports only full_pool + per-group aggregated outputs.
  - No FDR/gene-list outputs in randomize phase.
"""

import argparse
import os
import time
import logging
from tqdm import tqdm

from src.aggregation.concatenate import aggregate as aggregate_fn
from src.randomization.randomize import main as randomize_fn


def timed_execution(func, args, description):
    print(f"\n=== {description.upper()} ===")
    start_time = time.time()
    with tqdm(total=1, desc=description, bar_format="{l_bar}{bar} [elapsed: {elapsed}]") as pbar:
        func(args)
        pbar.update(1)
    print(f"{description} completed in {time.time() - start_time:.2f} seconds.\n")


def main():
    parser = argparse.ArgumentParser(
        description="CT_ACCUMULATION: CAAS gene accumulation randomization pipeline",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--tool", default="aggregate", choices=["aggregate", "randomize", "both"],
                        help="Phase(s) to run (default: aggregate)")
    parser.add_argument("--log-level", default="INFO",
                        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    parser.add_argument("--output-prefix", default="accumulation",
                        help="Prefix for output files (default: 'accumulation')")

    # Aggregation args
    parser.add_argument("--alignment-dir",    help="Directory containing alignment files")
    parser.add_argument("--genomic-info",     help="TSV: gene, chr, start, end, length")
    parser.add_argument("--species-list",     help="Traitfile (3-col, no header: species trait pair)")
    parser.add_argument("--metadata-caas",    help="Meta-CAAS file (global_meta_caas.tsv or original format)")
    parser.add_argument("--bg-caas",          help="Cleaned background gene list (no header)")
    parser.add_argument("--alignment-format", default="phylip-relaxed")

    # Randomization args
    parser.add_argument("--global-csv",                 help="Path to _global.csv from aggregation (contains masked + iscaas)")
    parser.add_argument("--caas-csv",                   help="Meta-CAAS file (global_meta_caas.tsv or CAAS CSV)")
    parser.add_argument("--randomization-type",         choices=["naive", "cons_decile"])
    parser.add_argument("--n-randomizations",           type=int, default=10000)
    parser.add_argument("--workers",                    type=int, default=None)
    parser.add_argument("--compress",                   action="store_true")
    parser.add_argument("--export-individual-rand",     action="store_true")
    parser.add_argument("--decile-bins",                type=str, default=None)
    parser.add_argument("--global-seed",                type=int, default=None)
    parser.add_argument("--precompute-masks",           dest="precompute_masks", action="store_true")
    parser.add_argument("--no-precompute-masks",        dest="precompute_masks", action="store_false")
    parser.set_defaults(precompute_masks=True)

    args = parser.parse_args()

    level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Starting CT_ACCUMULATION pipeline")

    print("Selected options:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    if args.tool in ["aggregate", "both"]:
        missing = [p for p in ["alignment_dir", "genomic_info", "species_list", "metadata_caas", "bg_caas"]
                   if getattr(args, p) is None]
        if missing:
            parser.error(f"Aggregation requires: {', '.join(missing)}")

    if args.tool in ["randomize", "both"]:
        missing = [p for p in ["global_csv", "caas_csv", "randomization_type"]
                   if getattr(args, p) is None]
        if missing:
            parser.error(f"Randomization requires: {', '.join(missing)}")

    try:
        if args.tool in ["aggregate", "both"]:
            agg_args = argparse.Namespace(
                alignment_dir=args.alignment_dir,
                alignment_format=args.alignment_format,
                genomic_info=args.genomic_info,
                species_list=args.species_list,
                metadata_caas=args.metadata_caas,
                bg_caas=args.bg_caas,
                output_prefix=args.output_prefix,
                log_level=args.log_level,
            )
            timed_execution(aggregate_fn, agg_args, "Aggregation Phase")

        if args.tool in ["randomize", "both"]:
            rand_args = argparse.Namespace(
                global_csv=args.global_csv,
                caas_csv=args.caas_csv,
                output_prefix=args.output_prefix,
                randomization_type=args.randomization_type,
                n_randomizations=args.n_randomizations,
                workers=args.workers,
                compress=args.compress,
                export_individual_rand=args.export_individual_rand,
                decile_bins=args.decile_bins,
                global_seed=args.global_seed,
                precompute_masks=args.precompute_masks,
                log_level=args.log_level,
            )
            timed_execution(randomize_fn, rand_args, "Randomization Phase")

    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

    logging.info("CT_ACCUMULATION pipeline completed successfully")


if __name__ == "__main__":
    main()
