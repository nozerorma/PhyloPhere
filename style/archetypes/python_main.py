#!/usr/bin/env python3
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: main.py
#

"""
CT_ACCUMULATION: Entry point for the CAAS gene-level accumulation pipeline.

Orchestrates two sequential phases:
  1. aggregate  — builds the global position matrix and per-group position CSVs
                  from filtered discovery output and alignment metadata.
  2. randomize  — permutation test estimating per-gene CAAS accumulation
                  significance against a random background.

Called by:  CT_ACCUMULATION Nextflow process (ctacc_run.nf → python main.py ...)
Usage:
    python main.py --tool [aggregate|randomize|both] [OPTIONS]

Run `python main.py --help` for the full option list.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import os
import time

# ── Third-party ───────────────────────────────────────────────────────────────
from tqdm import tqdm

# ── Package ───────────────────────────────────────────────────────────────────
from src.aggregation.concatenate import aggregate as aggregate_fn
from src.randomization.randomize import main as randomize_fn


# ── Helpers ───────────────────────────────────────────────────────────────────


def _timed(func, args, label: str) -> None:
    print(f"\n=== {label.upper()} ===", flush=True)
    t0 = time.time()
    with tqdm(total=1, desc=label, bar_format="{l_bar}{bar} [elapsed: {elapsed}]") as bar:
        func(args)
        bar.update(1)
    print(f"{label} completed in {time.time() - t0:.2f}s\n")


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CT_ACCUMULATION: CAAS gene accumulation randomization pipeline",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--tool", default="aggregate",
                   choices=["aggregate", "randomize", "both"],
                   help="Phase(s) to run (default: aggregate)")
    p.add_argument("--log-level", default="INFO",
                   choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    p.add_argument("--output-prefix", default="accumulation",
                   help="Prefix for all output files (default: 'accumulation')")

    # Aggregation phase
    agg = p.add_argument_group("aggregation")
    agg.add_argument("--alignment-dir",   help="Directory of alignment FASTA files")
    agg.add_argument("--genomic-info",    help="TSV: gene, chr, start, end, length")
    agg.add_argument("--species-list",    help="Traitfile (3-col, no header: species trait pair)")
    agg.add_argument("--metadata-caas",   help="Meta-CAAS file from CT_POSTPROC")

    # Randomization phase
    rnd = p.add_argument_group("randomization")
    rnd.add_argument("--n-permutations",  type=int, default=1000)
    rnd.add_argument("--change-side",     choices=["top", "bottom", "all"], default="top")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    if args.tool in ("aggregate", "both"):
        _timed(aggregate_fn, args, "Aggregation")
    if args.tool in ("randomize", "both"):
        _timed(randomize_fn, args, "Randomization")


if __name__ == "__main__":
    main()
