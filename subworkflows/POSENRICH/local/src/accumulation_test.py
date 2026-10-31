#!/usr/bin/env python3
"""Two explicit accumulation models, separately for top/bottom/all: uniform opportunity and canonical cycles.

The fixed family is every cleaned-background gene with tested positions. Each
cycle contributes its entire gene-level count, preserving dependence between
positions and the upstream cycle gene filter. The empirical test is unconditional over cycle counts. The hypergeometric test
conditions on the observed total K and assumes uniformly exchangeable positions
in fixed U. Neither test establishes positive selection.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import scipy.sparse as sp
from posenrich_contract import read_table, verify_universe, finish_manifest, sha256_file
from posenrich_within_gene import load_universe, read_block, benjamini_hochberg


def run(args):
    m = verify_universe(args.universe, args.universe_manifest, args.blocks)
    table = read_table(args.universe)
    U = load_universe(args.universe)
    genes = m['genes']
    gid = {g: i for i, g in enumerate(genes)}
    pos_genes = np.array([gid[g] for g in table.gene], dtype=np.int64)
    n_pool = np.bincount(pos_genes, minlength=len(genes))
    flag = {'top': 'observed_top', 'bottom': 'observed_bottom', 'both': 'observed_positive'}[args.change_side]
    actual = np.bincount(pos_genes, weights=table[flag].astype(int), minlength=len(genes)).astype(np.int64)
    # read_block uses U's sorted order, which need not be the background file order.
    ug = np.array([gid[U.genes[g]] for g in U.gene_of], dtype=np.int64)
    G = sp.csr_matrix((np.ones(len(ug)), (np.arange(len(ug)), ug)), shape=(len(U.keys), len(genes)))
    ge = np.zeros(len(genes), dtype=np.int64)
    total, total2, N = np.zeros(len(genes)), np.zeros(len(genes)), 0
    side = 'all' if args.change_side == 'both' else args.change_side
    for block in args.blocks:
        cycles, W = read_block(block, U, side=side)
        W.data = (W.data > 0).astype(float)
        W.eliminate_zeros()
        counts = (W @ G).toarray()
        ge += (counts >= actual).sum(axis=0)
        total += counts.sum(axis=0)
        total2 += (counts**2).sum(axis=0)
        N += len(cycles)
    if N != m['null']['N'] or N < 1:
        raise ValueError('Missing or incompatible null cycle denominator')
    eligible = n_pool > 0
    p = np.where(eligible, (ge + 1) / (N + 1), np.nan)
    mean = total / N
    sd = np.sqrt(np.maximum(total2 / N - mean**2, 0))
    M, K = int(n_pool.sum()), int(actual.sum())
    analytic = np.where(eligible, hypergeom.sf(actual - 1, M, n_pool, K), np.nan) if M else np.full(len(genes), np.nan)
    analytic_mean = K * n_pool / M if M else np.zeros(len(genes))
    analytic_sd = np.sqrt(K * (n_pool/M) * (1-n_pool/M) * (M-K)/(M-1)) if M>1 else np.zeros(len(genes))
    out = pd.DataFrame(dict(Gene=genes, ActualCount=actual, PoolPositions=n_pool,
                            MeanNull=analytic_mean, SD=analytic_sd, PValue=analytic, QValue=benjamini_hochberg(analytic),
                            EmpiricalMeanNull=mean, EmpiricalSD=sd, EmpiricalP=p, EmpiricalQ=benjamini_hochberg(p), N_cycles=N,
                            status=np.where(eligible, 'tested', 'no_tested_positions')))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL)
    finish_manifest(Path(args.output).parent, Path(args.output).name+'.manifest.json',
                    dict(component='accumulation', contract_version=2, direction=side,
                         inputs={'universe': sha256_file(args.universe),
                                 'universe_manifest': sha256_file(args.universe_manifest)},
                         N_cycles=N, family_size=int(eligible.sum()),
                         statistic='distinct positive positions per gene',
                         calibration={'hypergeometric':'upper tail Hyper(M, n_gene, K); uniform opportunity model',
                                      'empirical':'(1 + cycles with count >= observed)/(N + 1)'},
                         null_blocks=m['null']['blocks']), output_names=[Path(args.output).name])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--universe', required=True)
    p.add_argument('--universe-manifest', required=True)
    p.add_argument('--blocks', nargs='+', required=True)
    p.add_argument('--change-side', choices=['top', 'bottom', 'both'], default='both')
    p.add_argument('--output', required=True)
    run(p.parse_args())


if __name__ == '__main__':
    main()
