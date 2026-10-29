#!/usr/bin/env python3
# posenrich_within_gene.py — Level 2 of POSENRICH: concentration of signal inside genes, calibrated by the null cycles.
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""Does a term hold more positive positions than the rest of the positions of the same genes?

Estimand. For a term S and a gene g that has positions inside and outside S,
let n_g be the positions of g in U, m_g those in S and k_g the positive positions of g. If the positives of g were
spread evenly over its positions, k_g m_g / n_g of them would fall in S. The residual of the gene is

    r_g = (positives of g in S)  -  k_g m_g / n_g

and the statistic of the term is the sum of r_g over its eligible genes. Every gene is compared with itself, so
what the gene accumulates (level 1) does not enter the term. Observed and null cycles go through the same code
with the margins k_g of their own cycle, so detectability per position, co-occurrence of detections inside a gene
and the gene filter of the null are in the null.

Calibration: p = (1 + #{cycles with statistic >= observed}) / (N + 1), one-sided (concentration). Blocks of cycles
are read once and only the exceedance counts cross blocks, so the result does not depend on how cycles are split.

Eligibility is fixed by the annotation, never by the observed positives. A gene is eligible for a term when it has
members and non-members in U. A term is tested when it has members in at least --min-genes genes of U (default 2:
a term in one gene is a test inside that gene, not a concentration across genes); the rest are listed with the
reason and carry no p.

Hold-out mode (--holdout K): K null cycles play the observed one and the exceedances come from the other cycles.
The false-positive rate over all (term, cycle) pairs checks the procedure that is applied to the observed.

Descriptive effect sizes (no p): frequency, the proportion of positives among the positions of the term against the
same proportion in the rest of its eligible genes; intensity, the difference of mean score between the positives
inside and outside the term.

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
MIN_GENES = 2                        # genes of U a term needs members in to be tested
ALPHAS = (0.001, 0.01, 0.05)


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


def load_observed_scores(path, U):
    """Score per position of U: best-side positive CAAS score, 0 elsewhere (positive = score > 0)."""
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
    # the scores are only read for their sign: read_positives turns them into 0/1
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


def read_positives(block, U):
    """(cycle names, CSR matrix cycles x positions of U with 1 where the null cycle scores the position)."""
    cycles, W = read_block(block, U)
    W.data = (W.data > 0).astype(np.float64)
    W.eliminate_zeros()
    return cycles, W


# ── Residuals and statistic ───────────────────────────────────────────────────


def residual_chunks(W, design, pair_idx, chunk):
    """Yield (start, stop, R) with R[c, j] = r of pair pair_idx[start + j] in row c of W (cycles x pairs, dense)."""
    gene_positives = (W @ design.G).toarray()
    ratio = design.pair_m / design.pair_n
    for a in range(0, len(pair_idx), chunk):
        sel = pair_idx[a:a + chunk]
        A = (W @ design.M[:, sel]).toarray()
        yield a, a + len(sel), A - gene_positives[:, design.pair_gene[sel]] * ratio[sel]


def term_statistic(W, design, pair_idx, pair_to_term, chunk):
    """Statistic of every term for every row of W: array rows x terms (sum of r over the eligible genes)."""
    T = np.zeros((W.shape[0], len(design.terms)))
    for a, b, R in residual_chunks(W, design, pair_idx, chunk):
        T += R @ pair_to_term[a:b]
    return T


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


def eligibility(design, min_genes):
    """(testable term, eligible pair) masks. A term is tested when it has members in >= min_genes genes of U; a gene
    is eligible for it when it also has non-members. Both depend on the annotation and U only."""
    genes_per_term = np.bincount(design.pair_term, minlength=len(design.terms))
    term_ok = genes_per_term >= min_genes
    return term_ok, (design.pair_m < design.pair_n) & term_ok[design.pair_term]


def run(U, design, blocks, observed_score, min_genes=MIN_GENES, holdout=0, holdout_seed=1, chunk=20000):
    """Count the null cycles whose statistic reaches the target (the observed vector or the held-out cycles)."""
    cycle_lists = [block_cycles(b) for b in blocks]
    all_cycles = sorted(c for cs in cycle_lists for c in cs)
    if len(all_cycles) != len(set(all_cycles)):
        raise ValueError("a cycle appears in more than one block")
    held = set()
    if holdout:
        if holdout >= len(all_cycles):
            raise ValueError("--holdout must leave null cycles for the calibration")
        held = set(np.random.default_rng(holdout_seed).choice(all_cycles, size=holdout, replace=False).tolist())
    term_ok, eligible_pair = eligibility(design, min_genes)
    eligible = np.flatnonzero(eligible_pair)
    pair_to_term = sp.csr_matrix((np.ones(len(eligible)), (np.arange(len(eligible)), design.pair_term[eligible])),
                                 shape=(len(eligible), len(design.terms)))

    # targets: the observed vector, or the held-out cycles (read first: they must be known to count exceedances)
    if holdout:
        rows, target_names = [], []
        for block in blocks:
            cycles, W = read_positives(block, U)
            is_held = np.array([c in held for c in cycles])
            if is_held.any():
                rows.append(W[np.flatnonzero(is_held)])
                target_names += [c for c, h in zip(cycles, is_held) if h]
        Wt = sp.vstack(rows).tocsr()
    else:
        Wt = sp.csr_matrix((observed_score > 0).astype(np.float64).reshape(1, -1))
        target_names = ['observed']
    Tt = term_statistic(Wt, design, eligible, pair_to_term, chunk)

    ge = np.zeros(Tt.shape, dtype=np.int64)
    n_null = 0
    for block in blocks:
        cycles, W = read_positives(block, U)
        is_held = np.array([c in held for c in cycles])
        Wr = W[np.flatnonzero(~is_held)]
        n_null += Wr.shape[0]
        if Wr.shape[0] == 0 or len(eligible) == 0:
            continue
        ge += count_at_least(term_statistic(Wr, design, eligible, pair_to_term, chunk), Tt)
    return {'eligible': eligible, 'term_ok': term_ok, 'targets': target_names, 'T': Tt, 'ge': ge,
            'n_null': n_null, 'n_cycles_total': len(all_cycles), 'held': sorted(held)}


# ── Outputs ───────────────────────────────────────────────────────────────────


def results_tables(U, design, out, observed_score):
    eligible = out['eligible']
    n_terms = len(design.terms)
    eligible_pair = np.zeros(design.n_pairs, dtype=bool)
    eligible_pair[eligible] = True
    structural = design.pair_m < design.pair_n
    term_ok = out['term_ok']

    # observed quantities per pair
    positive = (observed_score > 0).astype(np.float64)
    pos_S = (sp.csr_matrix(positive.reshape(1, -1)) @ design.M).toarray().ravel()    # positives among the members
    pos_pair = (positive @ design.G)[design.pair_gene]                                # positives of the gene
    score_S = (sp.csr_matrix(observed_score.reshape(1, -1)) @ design.M).toarray().ravel()
    score_pair = (observed_score @ design.G)[design.pair_gene]
    expected = pos_pair * design.pair_m / design.pair_n

    reason = np.where(~structural, 'no_nonmembers', np.where(~term_ok[design.pair_term], 'term_below_min_genes', ''))
    by_gene = pd.DataFrame({
        'collection_id': design.terms['collection_id'].to_numpy()[design.pair_term],
        'term_id': design.terms['term_id'].to_numpy()[design.pair_term],
        'gene': np.asarray(U.genes, dtype=object)[design.pair_gene],
        'n_positions_gene': design.pair_n, 'n_members': design.pair_m,
        'positives_gene': pos_pair, 'positives_members': pos_S, 'expected_members': expected,
        'residual': np.where(structural, pos_S - expected, np.nan),
        'eligible': eligible_pair.astype(int), 'excluded_reason': reason})

    n_all = np.bincount(design.pair_term, minlength=n_terms)
    n_elig = np.bincount(design.pair_term[eligible_pair], minlength=n_terms)
    testable = n_elig > 0
    p = np.where(testable, (1.0 + out['ge'][0]) / (out['n_null'] + 1.0), np.nan)

    def per_term(values):
        return np.bincount(design.pair_term[eligible_pair], weights=values[eligible_pair], minlength=n_terms)
    members, n_pos = per_term(design.pair_m.astype(float)), per_term(design.pair_n.astype(float))
    ps, pg, exp_S = per_term(pos_S), per_term(pos_pair), per_term(expected)
    ss, sg = per_term(score_S), per_term(score_pair)
    with np.errstate(invalid='ignore', divide='ignore'):
        freq_in = ps / members
        freq_out = (pg - ps) / (n_pos - members)
        mean_in = ss / ps
        mean_out = (sg - ss) / (pg - ps)
    results = pd.DataFrame({
        'collection_id': design.terms['collection_id'], 'term_id': design.terms['term_id'],
        'n_genes_with_members_in_U': n_all, 'n_genes_eligible': n_elig, 'n_genes_excluded': n_all - n_elig,
        'n_members_eligible': members.astype(np.int64), 'positives_members_eligible': ps.astype(np.int64),
        'expected_positives_members': np.where(testable, exp_S, np.nan),
        'statistic': np.where(testable, out['T'][0], np.nan), 'p': p,
        'frequency_members': freq_in, 'frequency_rest': freq_out, 'frequency_difference': freq_in - freq_out,
        'mean_score_members': mean_in, 'mean_score_rest': mean_out, 'intensity_difference': mean_in - mean_out,
        'n_null_cycles': out['n_null']})
    results['q'] = np.nan
    for _, idx in results.groupby('collection_id').groups.items():
        results.loc[idx, 'q'] = benjamini_hochberg(results.loc[idx, 'p'].to_numpy())
    results['status'] = np.where(testable, 'tested', np.where(~term_ok, 'below_min_genes', 'not_identifiable'))
    return results, by_gene


def holdout_table(design, out):
    """False-positive rate over all (term, held-out cycle) pairs, by collection and alpha."""
    testable = np.zeros(len(design.terms), dtype=bool)
    testable[np.unique(design.pair_term[out['eligible']])] = True
    P = (1.0 + out['ge']) / (out['n_null'] + 1.0)
    rows = []
    for collection, idx in design.terms.groupby('collection_id').groups.items():
        cols = np.asarray(idx)[testable[np.asarray(idx)]]
        if len(cols) == 0:
            continue
        block = P[:, cols]
        for alpha in ALPHAS:
            rows.append({'collection_id': collection, 'alpha': alpha, 'n_terms_tested': len(cols),
                         'n_holdout_cycles': block.shape[0], 'n_pairs': block.size,
                         'fpr': float((block <= alpha).mean())})
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
    ap.add_argument('--min-genes', type=int, default=MIN_GENES, help='genes of U a term needs members in to be tested')
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
    observed_score = load_observed_scores(a.observed, U)
    logging.info(f"U: {len(U.keys)} positions; {len(design.terms)} terms, {design.n_pairs} (term, gene) pairs; "
                 f"{int((observed_score > 0).sum())} observed positives")
    out = run(U, design, a.blocks, observed_score, a.min_genes, a.holdout, a.holdout_seed, a.pair_chunk)
    logging.info(f"null cycles used: {out['n_null']} of {out['n_cycles_total']}; eligible pairs: {len(out['eligible'])}")

    os.makedirs(a.output_dir, exist_ok=True)
    results, by_gene = results_tables(U, design, out, observed_score)
    if a.holdout:
        write_tsv(holdout_table(design, out), os.path.join(a.output_dir, 'holdout_fpr.tsv'))
    else:
        write_tsv(results, os.path.join(a.output_dir, 'within_gene_results.tsv'))
        write_tsv(by_gene, os.path.join(a.output_dir, 'within_gene_by_gene.tsv'))
    manifest = {'component': 'within_gene', 'contract_version': 2,
                'inputs': {k: {'path': os.path.abspath(v), 'sha256': sha256_file(v)}
                           for k, v in (('universe', a.universe), ('observed', a.observed),
                                        ('collections', os.path.join(a.registry_dir, 'collections.tsv')))},
                'gene_list_check': gene_list_check(a.universe, a.registry_dir),
                'parameters': {'collections': collections, 'min_genes': a.min_genes,
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
