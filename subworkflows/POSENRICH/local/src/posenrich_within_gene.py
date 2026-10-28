#!/usr/bin/env python3
# posenrich_within_gene.py — Level 2 of POSENRICH: concentration of signal inside genes, calibrated by the null cycles.
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""Does a term hold more signal than the rest of the positions of the same genes?

Estimand (docs/POSENRICH_DISENO.md, section 6). Every position i of the universe U carries a weight w_i: its
best-side CAAS score when it is a positive detection, 0 otherwise. For a term S and a gene g that has positions
inside and outside S:

    r_g = sum of w_i over S ∩ g  -  W_g * m_g / n_g

with W_g the total weight of the gene, n_g its positions in U and m_g those in S. The statistic of the term is the
sum of r_g / s_g over the eligible genes (s_g: standard deviation of r_g over the null cycles); the plain sum of
r_g is a sensitivity. Observed and null cycles go through the same code with the margins of their own cycle, so
detectability per position, co-occurrence of detections inside a gene and the gene filter of the null are in the null.

Calibration: p = (1 + #{cycles with statistic >= observed}) / (N + 1), one-sided (concentration). A null cycle is
normalized with the s_g of the other cycles (leave-one-cycle-out, from the same counters), as the observed vector is,
because it does not enter its own s_g; with the s_g of all cycles a cycle would deflate its own statistic and the
observed one would look extreme. Blocks of cycles
are read twice and only additive counters cross blocks, so the result does not depend on how cycles are split:
  pass 1  per (term, gene): number of cycles, sum and sum of squares of r_g, cycles with r_g != 0   -> s_g
  pass 2  statistic of every cycle against the observed one                                         -> exceedances

Eligible genes: 0 < m_g < n_g (members and non-members in U) and a stable s_g (at least --min-nonzero-cycles null cycles
with r_g != 0 and s_g > 0). The rest are listed with their reason, never imputed.

Hold-out mode (--holdout K): K null cycles play the observed one, s_g and the exceedances come from the other
cycles. The false-positive rate over all (term, cycle) pairs checks the procedure that is applied to the observed.

Inputs:
  --universe       universe.tsv.gz of posenrich_universe.py (gene, position)
  --observed       position_scores.tsv or registry/position_evidence.tsv (Gene, Position, CAAS_score)
  --registry-dir   POSENRICH registry (positions.tsv, regions.tsv, sites.tsv, collections.tsv); its
                   informative_within_protein collections are the default ones
  --blocks         block directories of the null (positions.parquet: gene, cycle, position, all; cycles.txt)
Outputs (in --output_dir):
  within_gene_results.tsv  one row per (collection, term)
  within_gene_by_gene.tsv  one row per (term, gene) with members in U
  holdout_fpr.tsv          only with --holdout
  within_gene_manifest.json
"""

import argparse
import csv
import json
import logging
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import posenrich_members as pm  # noqa: E402
from posenrich_universe import POS_BITS, pair_keys, sha256_file  # noqa: E402

TIE_TOL = 1e-9                       # relative tolerance of ties between statistics
ZERO_TOL = 1e-12                     # |r| below it counts as zero
MIN_NONZERO_CYCLES = 30              # null cycles with r_g != 0 needed to estimate s_g (stability rule)
ALPHAS = (0.001, 0.01, 0.05)
STATS = ('normalized', 'simple')


# ── Inputs ────────────────────────────────────────────────────────────────────


@dataclass
class Universe:
    genes: list            # genes of U, index = gene id
    keys: np.ndarray       # sorted int64 keys (gene id, position) of U; the row index is the position index
    n_by_gene: np.ndarray  # positions of U per gene id

    @property
    def gene_of(self):
        return self.keys >> POS_BITS


def load_universe(path):
    df = pd.read_csv(path, sep='\t', usecols=['gene', 'position'], dtype={'gene': str, 'position': np.int64})
    genes = sorted(df['gene'].unique())
    ids = {g: i for i, g in enumerate(genes)}
    keys = np.unique(pair_keys(df['gene'].map(ids).to_numpy(), df['position'].to_numpy()))
    return Universe(genes, keys, np.bincount(keys >> POS_BITS, minlength=len(genes)))


def locate(U, genes, positions):
    """Index in U of each (gene, position); -1 when the position is not in U."""
    ids = pd.Series(genes).map({g: i for i, g in enumerate(U.genes)})
    ok = ids.notna().to_numpy()
    idx = np.full(len(ok), -1, dtype=np.int64)
    if ok.any():
        keys = pair_keys(ids[ok].to_numpy(dtype=np.int64), np.asarray(positions)[ok])
        at = np.minimum(np.searchsorted(U.keys, keys), len(U.keys) - 1)
        hit = U.keys[at] == keys
        idx[np.flatnonzero(ok)[hit]] = at[hit]
    return idx


def load_observed_weights(path, U):
    """Weight per position of U: best-side positive CAAS score, 0 elsewhere."""
    df = pd.read_csv(path, sep='\t', usecols=['Gene', 'Position', 'CAAS_score'], dtype=str, na_values=['NA'])
    df['score'] = pd.to_numeric(df['CAAS_score'], errors='coerce')
    df = df[np.isfinite(df['score']) & (df['score'] > 0)]
    best = df.assign(Position=df['Position'].astype(np.int64)).groupby(['Gene', 'Position'], as_index=False)['score'].max()
    idx = locate(U, best['Gene'].tolist(), best['Position'].to_numpy())
    if (idx < 0).any():
        logging.warning(f"{int((idx < 0).sum())} observed positive positions are not in U and are ignored")
    w = np.zeros(len(U.keys))
    w[idx[idx >= 0]] = best['score'].to_numpy()[idx >= 0]
    return w


def default_collections(registry_dir):
    """Collections whose terms are informative within a protein (collections.tsv of the registry)."""
    df = pd.read_csv(Path(registry_dir) / 'collections.tsv', sep='\t', dtype=str)
    return sorted(df.loc[df['informative_within_protein'] == '1', 'collection_id'].unique())


def load_members(registry_dir, collections, U):
    """DataFrame (collection_id, term_id, uidx): members of the terms of `collections` that lie in U.

    Collections are expanded one at a time from the compact registry (posenrich_members.py) and cut to U as they come.
    """
    positions, regions, sites = pm.read_tables(registry_dir, collections)
    parts = []
    for cid, frame in pm.iter_members(positions, regions, sites, collections):
        idx = locate(U, frame['Gene'].tolist(), frame['Position'].to_numpy())
        keep = idx >= 0
        parts.append(pd.DataFrame({'collection_id': cid, 'term_id': frame['term_id'].to_numpy()[keep],
                                   'uidx': idx[keep]}))
    if not parts:
        return pd.DataFrame({'collection_id': [], 'term_id': [], 'uidx': np.array([], dtype=np.int64)})
    return pd.concat(parts, ignore_index=True).drop_duplicates()


@dataclass
class Design:
    """Terms, (term, gene) pairs and the sparse maps between positions, pairs and terms."""
    terms: pd.DataFrame          # collection_id, term_id (row = term index)
    pair_term: np.ndarray        # term index per pair
    pair_gene: np.ndarray        # gene id per pair
    pair_m: np.ndarray           # positions of the term in the gene (in U)
    pair_n: np.ndarray           # positions of the gene in U
    M: sp.csc_matrix             # positions of U x pairs: 1 when the position is a member of the pair's term
    G: sp.csr_matrix             # positions of U x genes

    @property
    def n_pairs(self):
        return len(self.pair_term)


def build_design(members, U):
    terms = members[['collection_id', 'term_id']].drop_duplicates().sort_values(['collection_id', 'term_id'],
                                                                                kind='stable').reset_index(drop=True)
    term_index = pd.MultiIndex.from_frame(terms)
    t = term_index.get_indexer(pd.MultiIndex.from_frame(members[['collection_id', 'term_id']]))
    gene = U.gene_of[members['uidx'].to_numpy()]
    n_genes = len(U.genes)
    pkey = t.astype(np.int64) * n_genes + gene
    uniq, inverse, m = np.unique(pkey, return_inverse=True, return_counts=True)
    pair_term, pair_gene = uniq // n_genes, uniq % n_genes
    M = sp.csc_matrix((np.ones(len(inverse)), (members['uidx'].to_numpy(), inverse)), shape=(len(U.keys), len(uniq)))
    G = sp.csr_matrix((np.ones(len(U.keys)), (np.arange(len(U.keys)), U.gene_of)), shape=(len(U.keys), n_genes))
    return Design(terms, pair_term, pair_gene, m, U.n_by_gene[pair_gene], M, G)


# ── Null blocks ───────────────────────────────────────────────────────────────


def block_cycles(block):
    """Cycle names of a block: cycles.txt, else the distinct cycles of its positions table."""
    listing = Path(block) / 'cycles.txt'
    if listing.is_file():
        return [c for c in listing.read_text().splitlines() if c.strip()]
    return sorted(set(pq.read_table(Path(block) / 'positions.parquet', columns=['cycle'])['cycle'].to_pylist()))


def read_block(block, U):
    """(cycle names, CSR matrix cycles x positions of U with the null scores). One row per scored position and cycle."""
    cycles = block_cycles(block)
    table = pq.read_table(Path(block) / 'positions.parquet', columns=['gene', 'cycle', 'position', 'all'])
    gene = pc.index_in(table['gene'].combine_chunks(), value_set=pa.array(U.genes, pa.string()))
    cyc = pc.index_in(table['cycle'].combine_chunks(), value_set=pa.array(cycles, pa.string()))
    if pc.any(pc.is_null(cyc)).as_py():
        raise ValueError(f"{block}: positions of cycles that cycles.txt does not list")
    g = pc.fill_null(gene, -1).to_numpy().astype(np.int64)
    keep = g >= 0
    keys = pair_keys(g[keep], table['position'].to_numpy(zero_copy_only=False)[keep])
    at = np.minimum(np.searchsorted(U.keys, keys), len(U.keys) - 1)
    hit = U.keys[at] == keys
    rows = cyc.to_numpy().astype(np.int64)[keep][hit]
    cols = at[hit]
    score = table['all'].to_numpy(zero_copy_only=False)[keep][hit].astype(np.float64)
    W = sp.coo_matrix((score, (rows, cols)), shape=(len(cycles), len(U.keys))).tocsr()
    if W.nnz != len(score):
        raise ValueError(f"{block}: more than one row per (cycle, position)")
    return cycles, W


# ── Residuals and statistics ──────────────────────────────────────────────────


def residual_chunks(W, design, pair_idx, chunk):
    """Yield (start, stop, R) with R[c, j] = r of pair pair_idx[start + j] in row c of W (cycles x pairs, dense)."""
    gene_weight = (W @ design.G).toarray()
    ratio = design.pair_m / design.pair_n
    for a in range(0, len(pair_idx), chunk):
        sel = pair_idx[a:a + chunk]
        A = (W @ design.M[:, sel]).toarray()
        yield a, a + len(sel), A - gene_weight[:, design.pair_gene[sel]] * ratio[sel]


def loo_inverse_sd(R, s1, s2, n):
    """1 / s of every pair computed without the cycle of each row: R is cycles x pairs, s1/s2/n the counters of the
    pairs over all calibration cycles. Pairs whose leave-one-out variance is not positive scale to 0."""
    var = ((s2 - R * R) - (s1 - R) ** 2 / (n - 1)) / (n - 2)
    out = np.zeros_like(var)
    np.divide(1.0, np.sqrt(var), out=out, where=var > 0)
    return out


def term_statistics(W, design, pair_idx, inv_s, pair_to_term, chunk, loo=None):
    """(normalized, simple) statistics of every term for every row of W: arrays rows x terms.

    loo = (s1, s2, n) of the pairs of pair_idx, when the rows of W are calibration cycles (leave-one-cycle-out s);
    None for the observed vector or a held-out cycle, which use the s of all calibration cycles (inv_s).
    """
    n_terms = len(design.terms)
    Tn = np.zeros((W.shape[0], n_terms))
    Ts = np.zeros((W.shape[0], n_terms))
    for a, b, R in residual_chunks(W, design, pair_idx, chunk):
        B = pair_to_term[a:b]
        scale = inv_s[a:b] if loo is None else loo_inverse_sd(R, loo[0][a:b], loo[1][a:b], loo[2][a:b])
        Tn += (R * scale) @ B
        Ts += R @ B
    return Tn, Ts


def count_at_least(Tc, Tt):
    """Per term: number of rows of Tc with a statistic >= each target (ties within TIE_TOL count)."""
    S = np.sort(Tc, axis=0)
    thr = Tt - TIE_TOL * (1.0 + np.abs(Tt))
    out = np.empty(Tt.shape, dtype=np.int64)
    for j in range(Tc.shape[1]):
        out[:, j] = Tc.shape[0] - np.searchsorted(S[:, j], thr[:, j], side='left')
    return out


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        v = p[ok]
        order = np.argsort(v, kind='stable')
        ranked = v[order] * len(v) / (np.arange(len(v)) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(len(v))
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q


# ── Engine ────────────────────────────────────────────────────────────────────


@dataclass
class Fit:
    n_cycles: np.ndarray
    nonzero: np.ndarray
    sd: np.ndarray
    stable: np.ndarray           # structurally testable pair with a stable s_g


def run(U, design, blocks, observed_w, min_nonzero=MIN_NONZERO_CYCLES, holdout=0, holdout_seed=1, chunk=20000):
    """Run both passes. Returns a dict with the fit, the target statistics, the exceedance counts and N."""
    cycle_lists = [block_cycles(b) for b in blocks]
    all_cycles = sorted(c for cs in cycle_lists for c in cs)
    if len(all_cycles) != len(set(all_cycles)):
        raise ValueError("a cycle appears in more than one block")
    held = set()
    if holdout:
        if holdout >= len(all_cycles):
            raise ValueError("--holdout must leave null cycles for the calibration")
        held = set(np.random.default_rng(holdout_seed).choice(all_cycles, size=holdout, replace=False).tolist())
    structural = np.flatnonzero(design.pair_m < design.pair_n)

    # pass 1: counters of r over the calibration cycles; the held-out cycles are kept as targets
    n = np.zeros(design.n_pairs)
    s1 = np.zeros(design.n_pairs)
    s2 = np.zeros(design.n_pairs)
    nz = np.zeros(design.n_pairs)
    target_rows, target_names = [], []
    n_rest = 0
    for block in blocks:
        cycles, W = read_block(block, U)
        is_held = np.array([c in held for c in cycles])
        if is_held.any():
            target_rows.append(W[np.flatnonzero(is_held)])
            target_names += [c for c, h in zip(cycles, is_held) if h]
        Wr = W[np.flatnonzero(~is_held)]
        n_rest += Wr.shape[0]
        if Wr.shape[0] == 0 or len(structural) == 0:
            continue
        for a, b, R in residual_chunks(Wr, design, structural, chunk):
            sel = structural[a:b]
            n[sel] += R.shape[0]
            s1[sel] += R.sum(axis=0)
            s2[sel] += (R * R).sum(axis=0)
            nz[sel] += (np.abs(R) > ZERO_TOL).sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        var = (s2 - s1 * s1 / n) / (n - 1)
    sd = np.sqrt(np.clip(var, 0.0, None))
    stable = np.zeros(design.n_pairs, dtype=bool)
    stable[structural] = (nz[structural] >= min_nonzero) & (sd[structural] > 0)
    fit = Fit(n, nz, sd, stable)

    # targets: the observed vector, or the held-out cycles
    if holdout:
        Wt = sp.vstack(target_rows).tocsr()
    else:
        Wt = sp.csr_matrix(observed_w.reshape(1, -1))
        target_names = ['observed']
    eligible = np.flatnonzero(stable)
    pair_to_term = sp.csr_matrix((np.ones(len(eligible)), (np.arange(len(eligible)), design.pair_term[eligible])),
                                 shape=(len(eligible), len(design.terms)))
    inv_s = 1.0 / sd[eligible]
    loo = (s1[eligible], s2[eligible], n[eligible])
    Tn_t, Ts_t = term_statistics(Wt, design, eligible, inv_s, pair_to_term, chunk)

    # pass 2: exceedances of the calibration cycles over each target
    ge = {s: np.zeros(Tn_t.shape, dtype=np.int64) for s in STATS}
    for block in blocks:
        cycles, W = read_block(block, U)
        is_held = np.array([c in held for c in cycles])
        Wr = W[np.flatnonzero(~is_held)]
        if Wr.shape[0] == 0 or len(eligible) == 0:
            continue
        Tn, Ts = term_statistics(Wr, design, eligible, inv_s, pair_to_term, chunk, loo=loo)
        ge['normalized'] += count_at_least(Tn, Tn_t)
        ge['simple'] += count_at_least(Ts, Ts_t)
    return {'fit': fit, 'eligible': eligible, 'targets': target_names, 'T': {'normalized': Tn_t, 'simple': Ts_t},
            'ge': ge, 'n_null': n_rest, 'n_cycles_total': len(all_cycles), 'held': sorted(held)}


# ── Outputs ───────────────────────────────────────────────────────────────────


def term_p(counts, n_null, testable):
    p = (1.0 + counts) / (n_null + 1.0)
    return np.where(testable, p, np.nan)


def results_tables(U, design, out, observed_w):
    fit, eligible = out['fit'], out['eligible']
    n_terms = len(design.terms)
    eligible_pair = np.zeros(design.n_pairs, dtype=bool)
    eligible_pair[eligible] = True
    structural = design.pair_m < design.pair_n

    # observed quantities per pair
    wrow = sp.csr_matrix(observed_w.reshape(1, -1))
    prow = sp.csr_matrix((observed_w > 0).astype(float).reshape(1, -1))
    A = (wrow @ design.M).toarray().ravel()                      # weight of the members
    pos_S = (prow @ design.M).toarray().ravel()                  # positives among the members
    W_pair = (wrow @ design.G).toarray().ravel()[design.pair_gene]   # weight of the gene
    pos_pair = (prow @ design.G).toarray().ravel()[design.pair_gene]  # positives of the gene
    r_obs = A - W_pair * design.pair_m / design.pair_n

    reason = np.where(~structural, 'no_nonmembers',
                      np.where(~fit.stable, 'unstable_sd', ''))
    by_gene = pd.DataFrame({
        'collection_id': design.terms['collection_id'].to_numpy()[design.pair_term],
        'term_id': design.terms['term_id'].to_numpy()[design.pair_term],
        'gene': np.asarray(U.genes, dtype=object)[design.pair_gene],
        'n_positions_gene': design.pair_n, 'n_members': design.pair_m,
        'positives_gene': pos_pair, 'positives_members': pos_S,
        'weight_gene': W_pair, 'weight_members': A,
        'residual': np.where(structural, r_obs, np.nan),
        'sd_null': np.where(structural, fit.sd, np.nan),
        'cycles_nonzero': fit.nonzero.astype(np.int64),
        'eligible': eligible_pair.astype(int), 'excluded_reason': reason})

    p = {}
    testable = np.zeros(n_terms, dtype=bool)
    testable[np.unique(design.pair_term[eligible])] = True
    for stat in STATS:
        p[stat] = term_p(out['ge'][stat][0], out['n_null'], testable)

    et = design.pair_term[eligible_pair]
    def per_term(values):
        return np.bincount(design.pair_term[eligible_pair], weights=values[eligible_pair], minlength=n_terms)
    n_elig = np.bincount(et, minlength=n_terms)
    n_all = np.bincount(design.pair_term, minlength=n_terms)
    members = per_term(design.pair_m.astype(float))
    n_pos = per_term(design.pair_n.astype(float))
    ps, pg = per_term(pos_S), per_term(pos_pair)
    ws, wg = per_term(A), per_term(W_pair)
    with np.errstate(invalid='ignore', divide='ignore'):
        freq_in = ps / members
        freq_out = (pg - ps) / (n_pos - members)
        mean_in = ws / ps
        mean_out = (wg - ws) / (pg - ps)
    results = pd.DataFrame({
        'collection_id': design.terms['collection_id'], 'term_id': design.terms['term_id'],
        'n_genes_with_members_in_U': n_all, 'n_genes_eligible': n_elig, 'n_genes_excluded': n_all - n_elig,
        'n_members_eligible': members.astype(np.int64), 'positives_members_eligible': ps.astype(np.int64),
        'stat_normalized': np.where(testable, out['T']['normalized'][0], np.nan),
        'stat_simple': np.where(testable, out['T']['simple'][0], np.nan),
        'p_normalized': p['normalized'], 'p_simple': p['simple'],
        'frequency_members': freq_in, 'frequency_rest': freq_out, 'frequency_difference': freq_in - freq_out,
        'mean_score_members': mean_in, 'mean_score_rest': mean_out, 'intensity_difference': mean_in - mean_out,
        'n_null_cycles': out['n_null']})
    for stat in STATS:
        results[f'q_{stat}'] = np.nan
        for _, idx in results.groupby('collection_id').groups.items():
            results.loc[idx, f'q_{stat}'] = benjamini_hochberg(results.loc[idx, f'p_{stat}'].to_numpy())
    results['status'] = np.where(testable, 'tested', 'not_identifiable')
    return results, by_gene


def holdout_table(design, out):
    """False-positive rate over all (term, held-out cycle) pairs, by collection, statistic and alpha."""
    fit, eligible = out['fit'], out['eligible']
    testable = np.zeros(len(design.terms), dtype=bool)
    testable[np.unique(design.pair_term[eligible])] = True
    rows = []
    for stat in STATS:
        P = (1.0 + out['ge'][stat]) / (out['n_null'] + 1.0)
        for collection, idx in design.terms.groupby('collection_id').groups.items():
            cols = np.asarray(idx)[testable[np.asarray(idx)]]
            if len(cols) == 0:
                continue
            block = P[:, cols]
            for alpha in ALPHAS:
                rows.append({'collection_id': collection, 'statistic': stat, 'alpha': alpha,
                             'n_terms_tested': len(cols), 'n_holdout_cycles': block.shape[0],
                             'n_pairs': block.size, 'fpr': float((block <= alpha).mean())})
    return pd.DataFrame(rows)


def gene_list_check(universe_path, registry_dir):
    """Compare the gene-list digest of the universe with the one the registry recorded (digest of the list as a set).

    status: match, mismatch, or not_recorded when either manifest is missing or lacks the digest.
    """
    def digest(path, *keys):
        try:
            with open(path) as f:
                node = json.load(f)
            for key in keys:
                node = node[key]
            return node
        except (OSError, KeyError, TypeError, ValueError):
            return None

    universe = digest(os.path.join(os.path.dirname(os.path.abspath(universe_path)), 'universe_manifest.json'), 'gene_set_sha256')
    registry = digest(os.path.join(registry_dir, 'manifest.json'), 'universe', 'gene_list_check', 'registry_digest')
    status = 'not_recorded' if not universe or not registry else ('match' if universe == registry else 'mismatch')
    if status == 'mismatch':
        logging.warning("the gene list of the universe differs from the one of the registry: the terms and the universe "
                        "do not describe the same genes")
    return {'status': status, 'universe_digest': universe, 'registry_digest': registry}


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


# ── CLI ───────────────────────────────────────────────────────────────────────


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--universe', required=True)
    ap.add_argument('--observed', required=True)
    ap.add_argument('--registry-dir', required=True)
    ap.add_argument('--collections', default=None, help='comma-separated collection ids (default: informative_within_protein = 1)')
    ap.add_argument('--blocks', nargs='+', required=True)
    ap.add_argument('--output-dir', required=True)
    ap.add_argument('--min-nonzero-cycles', type=int, default=MIN_NONZERO_CYCLES)
    ap.add_argument('--holdout', type=int, default=0)
    ap.add_argument('--holdout-seed', type=int, default=1)
    ap.add_argument('--pair-chunk', type=int, default=20000)
    ap.add_argument('--log-level', default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'])
    a = ap.parse_args(argv)
    logging.basicConfig(level=getattr(logging, a.log_level), format='%(asctime)s - %(levelname)s - %(message)s')
    t0 = time.time()

    U = load_universe(a.universe)
    collections = [c for c in a.collections.split(',') if c] if a.collections else default_collections(a.registry_dir)
    members = load_members(a.registry_dir, collections, U)
    if members.empty:
        raise ValueError(f"no member of the collections {collections} lies in U")
    design = build_design(members, U)
    observed_w = load_observed_weights(a.observed, U)
    logging.info(f"U: {len(U.keys)} positions; {len(design.terms)} terms, {design.n_pairs} (term, gene) pairs; "
                 f"{int((observed_w > 0).sum())} observed positives")
    out = run(U, design, a.blocks, observed_w, a.min_nonzero_cycles, a.holdout, a.holdout_seed, a.pair_chunk)
    logging.info(f"null cycles used: {out['n_null']} of {out['n_cycles_total']}; eligible pairs: {len(out['eligible'])}")

    os.makedirs(a.output_dir, exist_ok=True)
    results, by_gene = results_tables(U, design, out, observed_w)
    if a.holdout:
        write_tsv(holdout_table(design, out), os.path.join(a.output_dir, 'holdout_fpr.tsv'))
    else:
        write_tsv(results, os.path.join(a.output_dir, 'within_gene_results.tsv'))
        write_tsv(by_gene, os.path.join(a.output_dir, 'within_gene_by_gene.tsv'))
    manifest = {'component': 'within_gene', 'contract_version': 1,
                'inputs': {k: {'path': os.path.abspath(v), 'sha256': sha256_file(v)}
                           for k, v in (('universe', a.universe), ('observed', a.observed),
                                        ('collections', os.path.join(a.registry_dir, 'collections.tsv')))},
                'gene_list_check': gene_list_check(a.universe, a.registry_dir),
                'parameters': {'collections': collections, 'min_nonzero_cycles': a.min_nonzero_cycles,
                               'holdout': a.holdout, 'holdout_seed': a.holdout_seed, 'tie_tol': TIE_TOL},
                'counts': {'n_universe': int(len(U.keys)), 'n_terms': int(len(design.terms)),
                           'n_pairs': int(design.n_pairs), 'n_pairs_eligible': int(len(out['eligible'])),
                           'n_null_cycles': int(out['n_null']), 'n_blocks': len(a.blocks),
                           'n_held_out': len(out['held'])},
                'seconds': round(time.time() - t0, 1)}
    with open(os.path.join(a.output_dir, 'within_gene_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
