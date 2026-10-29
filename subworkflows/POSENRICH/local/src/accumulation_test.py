#!/usr/bin/env python3
# accumulation_test.py — Level 1 of POSENRICH: does a gene hold a larger share of the detections than in the null?
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""
Accumulation test, calibrated with the null cycles.

Statistic
---------
Let K be the number of observed positive positions of the universe U and c_g those of gene g. The share of the gene
is c_g / K. Each null cycle has its own number of positive positions K_b and its own count c_g(b) per gene, so the
same share c_g(b) / K_b is computed in every cycle:

    p_g = (1 + #{cycles with share_g(b) >= c_g / K}) / (N + 1)

The null carries what a hypergeometric test ignores: positions are not equally detectable, detections co-occur
inside a gene and genes differ in how often they are scored. Dividing by K_b removes the difference between the
total number of detections of the observed data and that of a cycle, which would otherwise tilt every gene.
Only genes with at least one observed positive are tested; BH runs over them. The floor of p is 1 / (N + 1).

The null stores one score per position and cycle (best side), so the test applies to detections on either side; the
observed counts on the top and bottom sides are reported as descriptive columns.

Called by:  POSENRICH Nextflow process (accumulation_test.nf ACCUMULATION_TEST)
Inputs:     --universe (universe.tsv.gz), --blocks (block directories of the null: positions.parquet, cycles.txt)
Outputs:    accumulation_results.tsv (one row per gene of U), accumulation_manifest.json
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import posenrich_within_gene as wg  # noqa: E402

TIE_TOL = wg.TIE_TOL


# ── Observed counts ───────────────────────────────────────────────────────────


def observed_counts(universe_path, U):
    """Positions, observed positives (either side, top, bottom) per gene of U, in the gene order of U."""
    df = pd.read_csv(universe_path, sep='\t', usecols=['gene', 'position', 'observed_positive', 'observed_top',
                                                       'observed_bottom'],
                     dtype={'gene': str, 'position': np.int64})
    idx = wg.locate(U, df['gene'].tolist(), df['position'].to_numpy())
    gene = U.gene_of[idx]
    count = {k: np.bincount(gene, weights=df[col].to_numpy(), minlength=len(U.genes)).astype(np.int64)
             for k, col in (('all', 'observed_positive'), ('top', 'observed_top'), ('bottom', 'observed_bottom'))}
    return count


# ── Null ──────────────────────────────────────────────────────────────────────


def null_shares(blocks, U, observed_share):
    """Run the blocks once. Returns the exceedance count of every gene and the sums of its null share.

    Only additive counters cross blocks, so the result does not depend on how the cycles are split.
    """
    G = sp.csr_matrix((np.ones(len(U.keys)), (np.arange(len(U.keys)), U.gene_of)), shape=(len(U.keys), len(U.genes)))
    n_genes = len(U.genes)
    ge = np.zeros(n_genes, dtype=np.int64)
    s1, s2 = np.zeros(n_genes), np.zeros(n_genes)
    n_cycles, n_empty = 0, 0
    totals = []
    thr = observed_share - TIE_TOL * (1.0 + np.abs(observed_share))
    for block in blocks:
        _, W = wg.read_positives(block, U)
        counts = (W @ G).toarray()
        K = counts.sum(axis=1)
        keep = K > 0
        n_empty += int((~keep).sum())
        share = counts[keep] / K[keep, None]
        ge += (share >= thr).sum(axis=0)
        s1 += share.sum(axis=0)
        s2 += (share * share).sum(axis=0)
        n_cycles += share.shape[0]
        totals.append(K[keep])
    return ge, s1, s2, n_cycles, n_empty, (np.concatenate(totals) if totals else np.zeros(0))


def run(universe_path, blocks, output_dir):
    t0 = time.time()
    U = wg.load_universe(universe_path)
    if len(U.keys) == 0:
        raise ValueError(f"{universe_path} has no position: the universe U cannot be defined.")
    count = observed_counts(universe_path, U)
    K = int(count['all'].sum())
    if K == 0:
        raise ValueError("no observed positive position in U: the accumulation test has nothing to test")
    share = count['all'] / K
    logging.info(f"|U| = {len(U.keys)} positions in {len(U.genes)} genes; {K} observed positive positions")

    ge, s1, s2, n_cycles, n_empty, totals = null_shares(blocks, U, share)
    if n_cycles == 0:
        raise ValueError("the blocks hold no null cycle with a positive position")
    mean = s1 / n_cycles
    sd = np.sqrt(np.clip(s2 / n_cycles - mean * mean, 0.0, None) * n_cycles / max(n_cycles - 1, 1))
    tested = count['all'] > 0
    p = np.where(tested, (1.0 + ge) / (n_cycles + 1.0), np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = np.where(mean > 0, share / mean, np.nan)
    results = pd.DataFrame({
        'gene': U.genes, 'n_positions': U.n_by_gene, 'positives': count['all'],
        'positives_top': count['top'], 'positives_bottom': count['bottom'],
        'share': share, 'null_share_mean': mean, 'null_share_sd': sd, 'share_ratio': ratio,
        'p': p, 'q': wg.benjamini_hochberg(p), 'n_null_cycles': n_cycles})

    os.makedirs(output_dir, exist_ok=True)
    wg.write_tsv(results, os.path.join(output_dir, 'accumulation_results.tsv'))
    manifest = {'component': 'accumulation', 'contract_version': 2,
                'inputs': {'universe': {'path': os.path.abspath(universe_path), 'sha256': wg.sha256_file(universe_path)}},
                'counts': {'n_universe': int(len(U.keys)), 'n_genes': int(len(U.genes)),
                           'n_genes_tested': int(tested.sum()), 'n_blocks': len(blocks),
                           'n_null_cycles': int(n_cycles), 'n_null_cycles_without_positives': n_empty,
                           'observed_positives': K,
                           'null_positives_per_cycle': {
                               'mean': float(totals.mean()), 'sd': float(totals.std(ddof=1)) if len(totals) > 1 else 0.0,
                               'min': int(totals.min()), 'median': float(np.median(totals)), 'max': int(totals.max()),
                               'observed_rank_fraction': float((totals >= K).mean())}},
                'seconds': round(time.time() - t0, 1)}
    with open(os.path.join(output_dir, 'accumulation_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    logging.info(f"genes tested: {int(tested.sum())} of {len(U.genes)}; null cycles: {n_cycles}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def main(argv=None):
    parser = argparse.ArgumentParser(description="POSENRICH level 1, accumulation: share of detections per gene, "
                                                 "calibrated with the null cycles",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--universe", required=True, help="universe.tsv.gz of posenrich_universe.py")
    parser.add_argument("--blocks", nargs='+', required=True, help="block directories of the null")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-level", default="INFO", choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"])
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format='%(asctime)s - %(levelname)s - %(message)s')
    run(args.universe, args.blocks, args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
