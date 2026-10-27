#!/usr/bin/env python3
# accumulation_test.py — Per-gene exact test of CAAS accumulation.
# PhyloPhere | subworkflows/CT_ACCUMULATION/local/

"""
Accumulation test: does a gene hold more CAAS positions than expected from its share of the position universe?

Universe U
----------
U is the set of positions an observed CAAS could come from, the same universe SCORING and POSENRICH use
(built by posenrich_universe.py): the positions CAAStools tested (background.output) in the genes of the
gene list (--bg-caas), kept when they score > 0 in at least one cycle and side of the filtered permulation
null (null_positions.parquet: one row per position the null scores, with the number of cycles that score it),
plus every observed positive position (a top or bottom row of --caas-csv with a positive asr_path_score).
A position the null never scores cannot produce a CAAS under the permulation labelings, so it is not part
of the opportunity of a gene.

Test
----
Per direction (top, bottom, both) let M be the number of distinct positions of U, K the number of distinct
observed CAAS positions of the direction (a position detected on both sides counts once in 'both'), n_g the
number of positions of gene g in U and obs_g its observed CAAS positions. Under the null the K observed
positions are a uniform draw without replacement from U, so the count of gene g is hypergeometric:
    p_g = P(X >= obs_g),  X ~ Hypergeometric(M, n_g, K)
    MeanNull_g = K n_g / M,  SD_g = sqrt(K (n_g/M) (1 - n_g/M) (M - K) / (M - 1))
The gene length enters through n_g. Genes without observed CAAS have p = 1 and are not tested; the
Benjamini-Hochberg adjustment over the tested genes is done downstream in accum_gene_lists.nf and
scoring_compute.R. Positions are treated as independent draws, which understates the variance of genes with
many positions in U.

Called by:  CT_ACCUMULATION Nextflow process (ctacc_run.nf CT_ACCUMULATION_TEST)
Inputs:     --caas-csv (filtered_discovery.tsv: one row per gene, position and side), --background-positions
            (caastools background.output), --bg-caas (cleaned_background_main.txt), --null-positions
            (null_positions.parquet)
Outputs:    --output CSV with columns Gene, ActualCount, PoolPositions, MeanNull, SD, PValue
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import os
import sys
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

# The universe U is built by the module POSENRICH shares (one definition for every consumer).
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ENRICHMENT' / 'local' / 'src'))
from posenrich_universe import POS_BITS, build_universe, load_null_positive, load_observed_positive, load_universe_genes  # noqa: E402

COLUMNS = ['Gene', 'ActualCount', 'PoolPositions', 'MeanNull', 'SD', 'PValue']


# ── Test ──────────────────────────────────────────────────────────────────────


def hypergeometric_test(obs, n_pool, M, K):
    """Per-gene p-value, null mean and null SD of the hypergeometric count.

    obs and n_pool are per-gene arrays (observed CAAS positions and positions of the gene in U), M the
    size of U and K the number of observed CAAS positions. Genes with obs = 0 get p = 1.
    """
    obs, n_pool = np.asarray(obs, dtype=np.int64), np.asarray(n_pool, dtype=np.int64)
    p = np.ones(len(obs))
    tested = obs > 0
    p[tested] = hypergeom.sf(obs[tested] - 1, M, n_pool[tested], K)
    share = n_pool / M
    mean = K * share
    sd = np.sqrt(K * share * (1 - share) * (M - K) / (M - 1)) if M > 1 else np.zeros(len(obs))
    return p, mean, sd


def run(args):
    """Run the exact accumulation test of one direction and write the results CSV.

    Reads from args: caas_csv, background_positions, bg_caas, null_positions, change_side, output.
    """
    genes = load_universe_genes(args.bg_caas)
    logging.info(f"Universe genes: {len(genes)}")

    null_positive = load_null_positive(args.null_positions)
    if null_positive.empty:
        raise ValueError(f"{args.null_positions} has no position the null scores (a null without cycles, for "
                         f"example --caas_full_perms 0): the universe U cannot be defined.")
    universe = build_universe(genes, null_positive, load_observed_positive(args.caas_csv), args.background_positions)
    m = universe.metrics
    logging.info(f"Null-positive in U: {m['n_null_positive_in_universe']}; |U| = {m['n_universe']}; tested set: "
                 f"{m['tested_set']}; {m.get('n_observed_outside_tested', 0)} observed positions lie outside the tested "
                 f"set and stay in U")
    M = int(universe.keys.size)
    n_genes = len(genes)
    n_pool = np.bincount(universe.keys >> POS_BITS, minlength=n_genes)

    # Observed positive positions of the direction: distinct positions, so a position detected on both sides
    # is one observed position in 'both'
    observed = {'top': universe.top, 'bottom': universe.bottom, 'both': universe.observed}[args.change_side]
    K = int(observed.size)
    actual = np.bincount(observed >> POS_BITS, minlength=n_genes)
    logging.info(f"Direction '{args.change_side}': {K} distinct observed positive positions; "
                 f"fraction of U observed = {K / M:.4f}")

    p, mean, sd = hypergeometric_test(actual, n_pool, M, K)
    logging.info(f"Genes tested (observed positives > 0): {int((actual > 0).sum())} of {n_genes}")

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    pd.DataFrame({'Gene': genes, 'ActualCount': actual, 'PoolPositions': n_pool,
                  'MeanNull': mean, 'SD': sd, 'PValue': p}, columns=COLUMNS).to_csv(args.output, index=False)
    logging.info(f"Results: {args.output}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="CT_ACCUMULATION: exact per-gene accumulation test",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--caas-csv", required=True, help="filtered_discovery.tsv (one row per gene, position and side)")
    parser.add_argument("--background-positions", required=True,
                        help="caastools background.output (gene<TAB>tested msa positions)")
    parser.add_argument("--bg-caas", required=True, help="cleaned_background_main.txt (genes that survive post-processing)")
    parser.add_argument("--null-positions", required=True,
                        help="null_positions.parquet from CAAS_CORE_MERGE_REDUCE, the filtered null that SCORING also uses")
    parser.add_argument("--change-side", default="both", choices=["top", "bottom", "both"],
                        help="Phenotype direction tested (default: both, each position counted once)")
    parser.add_argument("--output", required=True, help="Results CSV")
    parser.add_argument("--log-level", default="INFO", choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format='%(asctime)s - %(levelname)s - %(message)s')
    run(args)


if __name__ == "__main__":
    main()
