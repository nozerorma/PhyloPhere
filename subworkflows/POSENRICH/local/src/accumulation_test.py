#!/usr/bin/env python3
# accumulation_test.py — Per-gene exact test of CAAS accumulation.
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""
Accumulation test: does a gene hold more CAAS positions than expected from its share of the position universe?

Universe U
----------
U is read from universe.tsv.gz, the table posenrich_universe.py writes and level 2 also uses: one row per position
of U (gene, position, nd, observed_positive, observed_top, observed_bottom). Accumulation builds no universe of its
own.

Test
----
Per direction (top, bottom, both) let M be the number of rows of U, K the number of observed positive positions of
the direction (a position detected on both sides counts once in 'both'), n_g the number of positions of gene g in
U and obs_g its observed positive positions. Under the null the K observed positions are a uniform draw without
replacement from U, so the count of gene g is hypergeometric:
    p_g = P(X >= obs_g),  X ~ Hypergeometric(M, n_g, K)
    MeanNull_g = K n_g / M,  SD_g = sqrt(K (n_g/M) (1 - n_g/M) (M - K) / (M - 1))
The gene length enters through n_g. Genes of the gene list without a position in U (n_g = 0) have no row in the
output and are not tested. Genes without observed positives have p = 1 and are not tested either; the
Benjamini-Hochberg adjustment over the tested genes is done by the report. Positions are treated as equivalent
draws, which understates the variance of genes with many positions in U.

Called by:  POSENRICH Nextflow process (accumulation_test.nf ACCUMULATION_TEST)
Inputs:     --universe (universe.tsv.gz), --change-side (top, bottom or both)
Outputs:    --output CSV with columns Gene, ActualCount, PoolPositions, MeanNull, SD, PValue
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import os

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

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

    Reads from args: universe, change_side, output.
    """
    universe = pd.read_csv(args.universe, sep='\t', usecols=['gene', 'observed_positive', 'observed_top', 'observed_bottom'])
    M = len(universe)
    if M == 0:
        raise ValueError(f"{args.universe} has no position: the universe U cannot be defined.")
    # Genes in order of first appearance (the universe table is sorted by gene list order and position)
    genes, first = np.unique(universe['gene'].to_numpy(), return_index=True)
    genes = genes[np.argsort(first)]
    gene_id = pd.Categorical(universe['gene'], categories=genes).codes
    n_pool = np.bincount(gene_id, minlength=len(genes))
    logging.info(f"|U| = {M} positions in {len(genes)} genes")

    # Observed positive positions of the direction: distinct rows of U, so a position detected on both sides
    # is one observed position in 'both'
    flag = {'top': 'observed_top', 'bottom': 'observed_bottom', 'both': 'observed_positive'}[args.change_side]
    observed = universe[flag].to_numpy().astype(bool)
    K = int(observed.sum())
    actual = np.bincount(gene_id[observed], minlength=len(genes))
    logging.info(f"Direction '{args.change_side}': {K} distinct observed positive positions; "
                 f"fraction of U observed = {K / M:.4f}")

    p, mean, sd = hypergeometric_test(actual, n_pool, M, K)
    logging.info(f"Genes tested (observed positives > 0): {int((actual > 0).sum())} of {len(genes)}")

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    pd.DataFrame({'Gene': genes, 'ActualCount': actual, 'PoolPositions': n_pool,
                  'MeanNull': mean, 'SD': sd, 'PValue': p}, columns=COLUMNS).to_csv(args.output, index=False)
    logging.info(f"Results: {args.output}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="POSENRICH level 1, accumulation: exact per-gene accumulation test",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--universe", required=True, help="universe.tsv.gz of posenrich_universe.py")
    parser.add_argument("--change-side", default="both", choices=["top", "bottom", "both"],
                        help="Phenotype direction tested (default: both, each position counted once)")
    parser.add_argument("--output", required=True, help="Results CSV")
    parser.add_argument("--log-level", default="INFO", choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format='%(asctime)s - %(levelname)s - %(message)s')
    run(args)


if __name__ == "__main__":
    main()
