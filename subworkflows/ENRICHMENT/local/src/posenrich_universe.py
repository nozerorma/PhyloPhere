#!/usr/bin/env python3
# posenrich_universe.py — Universe U of positions shared by Accumulation and POSENRICH.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Build the universe U: the positions an observed CAAS could come from.

U = (positions CAAStools tested  ∩  positions that score > 0 in at least one cycle of the filtered null)
    restricted to the genes of the gene list, plus every observed positive position of those genes.

A position the null never scores cannot produce a CAAS under the permulation labelings, so it is not part of
the opportunity of a gene. An observed position the null never scores stays in U (its nd is 0).

The gene list is an explicit input and its hash goes to the manifest, so every consumer can check that it
used the same universe. When fewer than RATE_MIN of the null-positive positions lie in the tested set (the tested
file is in another coordinate space), the tested set is ignored, the same rule as SCORING (null_position_stats.py).

Observed positions are the rows with a positive score on the top or bottom side: `asr_path_score` of
filtered_discovery.tsv or `CAAS_score` of position_scores.tsv.

Used as a library by CT_ACCUMULATION (accumulation_test.py). As a command it writes:
  universe.tsv.gz        gene, position, nd (cycles that score the position), observed_positive (0/1)
  universe_manifest.json input hashes, counts and the tested-set decision
"""

import argparse
import hashlib
import json
import logging
import math
import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

POS_BITS = 32                  # position keys: gene id << 32 | position
RATE_MIN = 0.5                 # share of null-positive positions that must lie in the tested set to use it
SIDES = ('top', 'bottom')


def pair_keys(gene_ids, positions):
    """One int64 key per (gene id, position)."""
    return (np.asarray(gene_ids, dtype=np.int64) << POS_BITS) + np.asarray(positions, dtype=np.int64)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# ── Inputs ────────────────────────────────────────────────────────────────────


def load_universe_genes(path):
    """Genes of the gene list (one gene per line, optional 'Gene' header), in file order without repeats."""
    with open(path) as f:
        genes = [line.strip() for line in f if line.strip() and line.strip() != 'Gene']
    return list(dict.fromkeys(genes))


def read_tested(background_file, genes):
    """{gene: sorted unique int64 positions CAAStools tested}, for the genes in `genes`.

    background.output is one line per gene, `gene<TAB>comma-separated positions`, in the coordinate space of
    filtered_discovery's `Position`.
    """
    tested = {}
    with open(background_file) as f:
        for line in f:
            parts = line.rstrip('\n').rstrip('\r').split('\t')
            if len(parts) < 2 or parts[0] not in genes or parts[1] in ('', 'NULL', 'Position'):
                continue
            values = [int(p) for p in parts[1].split(',') if p.strip().lstrip('-').isdigit()]
            tested[parts[0]] = np.unique(np.array(values, dtype=np.int64))
    return tested


def load_null_positive(null_path):
    """DataFrame (gene, position, nd) of the positions the filtered null scores (nd > 0), from null_positions.parquet."""
    df = pd.read_parquet(null_path)
    missing = [c for c in ('gene', 'position', 'nd') if c not in df.columns]
    if missing:
        raise ValueError(f"{null_path} lacks column(s) {missing}: regenerate the CAAS permulation null (CAAS_CORE_MERGE_REDUCE).")
    df = df[['gene', 'position', 'nd']]
    df = df[df['nd'] > 0]
    return df.astype({'gene': str, 'position': np.int64, 'nd': np.int64}).reset_index(drop=True)


def is_positive_score(x):
    """True for a finite score > 0 (the rule of core.scores.position_score)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(v) and v > 0.0


def load_observed_positive(path):
    """DataFrame (gene, position, side) of the distinct observed positive rows of the top and bottom sides.

    `path` is filtered_discovery.tsv (score column asr_path_score) or position_scores.tsv (CAAS_score).
    """
    empty = pd.DataFrame({'gene': pd.Series(dtype=str), 'position': pd.Series(dtype=np.int64),
                          'side': pd.Series(dtype=str)})
    if os.path.getsize(path) == 0:
        logging.warning("observed positions file is empty: no positions passed the filter")
        return empty
    header = list(pd.read_csv(path, sep='\t', nrows=0).columns)
    score_col = next((c for c in ('asr_path_score', 'CAAS_score') if c in header), None)
    missing = [c for c in ('Gene', 'Position', 'side') if c not in header]
    if score_col is None or missing:
        raise ValueError(f"{path} must have Gene, Position, side and asr_path_score or CAAS_score "
                         f"(filtered_discovery.tsv or position_scores.tsv).")
    df = pd.read_csv(path, sep='\t', usecols=['Gene', 'Position', 'side', score_col])
    df['side'] = df['side'].astype(str).str.strip().str.lower()
    df = df[df['side'].isin(SIDES)]
    n_rows = len(df)
    df = df[df[score_col].map(is_positive_score)]
    if len(df) < n_rows:
        logging.info(f"{n_rows - len(df)} of {n_rows} observed rows have no positive score and are not observed positives")
    out = pd.DataFrame({'gene': df['Gene'].astype(str), 'position': df['Position'].astype(np.int64),
                        'side': df['side']})
    return out.drop_duplicates().reset_index(drop=True)


# ── Universe ──────────────────────────────────────────────────────────────────


@dataclass
class Universe:
    genes: list                      # gene list, index = gene id
    gene_to_id: dict
    keys: np.ndarray                 # sorted int64 keys (gene id, position) of U
    nd: np.ndarray                   # cycles that score each position of U (0 for observed positions the null never scores)
    top: np.ndarray                  # keys of the observed positive positions on the top side
    bottom: np.ndarray               # keys of the observed positive positions on the bottom side
    metrics: dict = field(default_factory=dict)

    @property
    def observed(self):
        """Keys of the distinct observed positive positions (either side)."""
        return np.union1d(self.top, self.bottom)

    def table(self):
        """DataFrame (gene, position, nd, observed_positive) of U."""
        observed = np.isin(self.keys, self.observed)
        return pd.DataFrame({'gene': np.asarray(self.genes, dtype=object)[self.keys >> POS_BITS],
                             'position': self.keys & ((1 << POS_BITS) - 1),
                             'nd': self.nd, 'observed_positive': observed.astype(np.int8)})


def _membership(frame, tested):
    """Boolean per row of `frame` (gene, position): the position is tested; and per row: its gene has a tested line."""
    in_tested = np.zeros(len(frame), dtype=bool)
    covered = np.zeros(len(frame), dtype=bool)
    positions = frame['position'].to_numpy()
    for gene, idx in frame.groupby('gene', sort=False).indices.items():
        if gene in tested:
            covered[idx] = True
            in_tested[idx] = np.isin(positions[idx], tested[gene])
    return in_tested, covered


def build_universe(genes, null_df, observed_df, background_file=None, rate_min=RATE_MIN):
    """U from the gene list, the null-positive positions, the observed positives and the tested positions."""
    genes = list(dict.fromkeys(genes))
    gene_to_id = {g: i for i, g in enumerate(genes)}
    metrics = {'n_genes': len(genes), 'n_null_positive_all': int(len(null_df))}

    in_universe = null_df['gene'].isin(gene_to_id).to_numpy()
    keep = in_universe.copy()
    tested = None
    metrics['tested_set'] = 'not_given'
    if background_file is not None:
        wanted = set(genes) | set(null_df['gene'].unique())
        tested = read_tested(background_file, wanted)
        in_tested, covered = _membership(null_df, tested)
        n_cov = int(covered.sum())
        rate = float(in_tested[covered].sum() / n_cov) if n_cov > 0 else 0.0
        metrics['tested_rate'] = rate
        if rate >= rate_min:
            keep = in_universe & in_tested
            metrics['tested_set'] = 'used'
            metrics['n_tested_positions'] = int(sum(len(v) for g, v in tested.items() if g in gene_to_id))
        else:
            metrics['tested_set'] = 'ignored_incompatible'
            logging.warning(f"tested positions ignored: only {100 * rate:.1f}% of the null-positive positions lie in them")
    null_u = null_df[keep]
    null_keys = pair_keys(null_u['gene'].map(gene_to_id).to_numpy(), null_u['position'].to_numpy())
    metrics['n_null_positive_in_universe'] = int(len(null_u))

    outside = observed_df[~observed_df['gene'].isin(gene_to_id)]
    if len(outside):
        logging.warning(f"{len(outside)} observed rows belong to genes outside the gene list and are dropped")
    obs = observed_df[observed_df['gene'].isin(gene_to_id)]
    top = np.unique(pair_keys(obs.loc[obs['side'] == 'top', 'gene'].map(gene_to_id).to_numpy(),
                              obs.loc[obs['side'] == 'top', 'position'].to_numpy()))
    bottom = np.unique(pair_keys(obs.loc[obs['side'] == 'bottom', 'gene'].map(gene_to_id).to_numpy(),
                                 obs.loc[obs['side'] == 'bottom', 'position'].to_numpy()))
    observed = np.union1d(top, bottom)
    metrics['n_observed_positive'] = int(observed.size)
    if tested is not None and metrics['tested_set'] == 'used':
        obs_frame = pd.DataFrame({'gene': np.asarray(genes, dtype=object)[observed >> POS_BITS],
                                  'position': observed & ((1 << POS_BITS) - 1)})
        obs_in_tested, _ = _membership(obs_frame, tested)
        metrics['n_observed_outside_tested'] = int((~obs_in_tested).sum())

    # nd of every position of U comes from all null rows of the list's genes: an observed position that the
    # null scores but CAAStools did not test keeps its nd.
    in_genes = null_df[in_universe]
    gene_keys = pair_keys(in_genes['gene'].map(gene_to_id).to_numpy(), in_genes['position'].to_numpy())
    order = np.argsort(gene_keys, kind='stable')
    sorted_keys = gene_keys[order]
    sorted_nd = in_genes['nd'].to_numpy()[order]
    keys = np.union1d(null_keys, observed)
    nd = np.zeros(len(keys), dtype=np.int64)
    at = np.minimum(np.searchsorted(sorted_keys, keys), max(len(sorted_keys) - 1, 0))
    found = (sorted_keys[at] == keys) if len(sorted_keys) else np.zeros(len(keys), dtype=bool)
    nd[found] = sorted_nd[at[found]]
    metrics['n_universe'] = int(keys.size)
    metrics['n_observed_not_scored_by_null'] = int((~np.isin(observed, sorted_keys)).sum())
    return Universe(genes, gene_to_id, keys, nd, top, bottom, metrics)


def manifest(universe, inputs):
    """Manifest dict: SHA-256 of every input file, counts and the tested-set decision."""
    return {'component': 'universe',
            'contract_version': 1,
            'inputs': {k: {'path': os.path.abspath(v), 'sha256': sha256_file(v)} for k, v in inputs.items() if v},
            'gene_list_sha256': sha256_file(inputs['background_genes']),
            'metrics': universe.metrics}


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--background-genes', required=True, help='gene list, one gene per line')
    p.add_argument('--background-positions', default=None, help='caastools background.output (gene<TAB>tested positions)')
    p.add_argument('--null-positions', required=True, help='null_positions.parquet')
    p.add_argument('--observed', required=True, help='filtered_discovery.tsv or position_scores.tsv')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--log-level', default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'])
    a = p.parse_args()
    logging.basicConfig(level=getattr(logging, a.log_level), format='%(asctime)s - %(levelname)s - %(message)s')

    genes = load_universe_genes(a.background_genes)
    null_df = load_null_positive(a.null_positions)
    if null_df.empty:
        raise ValueError(f"{a.null_positions} has no position the null scores: the universe U cannot be defined.")
    universe = build_universe(genes, null_df, load_observed_positive(a.observed), a.background_positions)

    os.makedirs(a.output_dir, exist_ok=True)
    universe.table().to_csv(os.path.join(a.output_dir, 'universe.tsv.gz'), sep='\t', index=False)
    inputs = {'background_genes': a.background_genes, 'background_positions': a.background_positions,
              'null_positions': a.null_positions, 'observed': a.observed}
    with open(os.path.join(a.output_dir, 'universe_manifest.json'), 'w') as f:
        json.dump(manifest(universe, inputs), f, indent=2, sort_keys=True)
        f.write('\n')
    logging.info(f"|U| = {universe.metrics['n_universe']} "
                 f"({universe.metrics['n_observed_positive']} observed positive positions)")


if __name__ == '__main__':
    main()
