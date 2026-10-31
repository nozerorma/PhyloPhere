#!/usr/bin/env python3
# posenrich_universe.py — Universe U of positions shared by Accumulation and POSENRICH.
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""Fixed positional opportunity: declared tested positions intersected with the cleaned gene background.

Observed and null detections annotate U; neither selects its members. Canonical null
blocks, their complete cycle roster and observed scientific provenance are verified
before any inferential consumer runs. Silent valid cycles remain in N.
"""

import argparse
import hashlib
import json
import logging
import math
import os
import sys
from pathlib import Path
from posenrich_contract import observed_table, identifier, position, block_inventory, finish_manifest, read_table
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

POS_BITS = 32                  # position keys: gene id << 32 | position
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


def gene_list_digest(genes):
    """SHA-256 of a gene list as a set: sorted unique names joined by newlines.

    SCORING records the same digest of its gene list (null_position_stats.py gene_list_identity), so the order of
    the file or a header line does not matter when two consumers compare lists.
    """
    return hashlib.sha256("\n".join(sorted(set(genes))).encode()).hexdigest()


# ── Inputs ────────────────────────────────────────────────────────────────────


def load_universe_genes(path):
    """Genes of the gene list (one gene per line, optional 'Gene' header), in file order without repeats."""
    genes = Path(path).read_text().splitlines()
    if genes and genes[0] in ('Gene', 'gene'):
        genes = genes[1:]
    for g in genes:
        identifier(g)
    if len(genes) != len(set(genes)):
        raise ValueError('Duplicate cleaned-background genes')
    return genes


def read_tested(background_file, genes):
    tested = {}
    for line in Path(background_file).read_text().splitlines():
        parts = line.split('\t')
        if len(parts) != 2 or identifier(parts[0]) in tested:
            raise ValueError('Invalid/duplicate tested gene')
        values = [] if parts[1] in ('', 'NULL') else [position(p) for p in parts[1].split(',')]
        if len(values) != len(set(values)):
            raise ValueError('Duplicate tested position')
        tested[parts[0]] = np.asarray(sorted(values), dtype=np.int64)
    return {g: v for g, v in tested.items() if g in genes}


def load_null_positive(null_path):
    df = pd.read_parquet(null_path)
    if not {'gene', 'position', 'nd'} <= set(df):
        raise ValueError('Missing aggregate null position columns')
    if df.duplicated(['gene', 'position']).any():
        raise ValueError('Duplicate aggregate null position')
    for g in df.gene:
        identifier(g)
    df['position'] = df['position'].map(position).astype(np.int64)
    if any(not str(n).isdigit() or int(n) < 1 for n in df.nd):
        raise ValueError('Invalid null detection counts')
    return df[['gene', 'position', 'nd']].astype({'nd': np.int64})


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
    df = observed_table(path)
    df = df[df.side.isin(SIDES)]
    return df[['Gene', 'Position', 'side']].rename(columns={'Gene': 'gene', 'Position': 'position'}).reset_index(drop=True)


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
        """DataFrame (gene, position, nd, observed_positive, observed_top, observed_bottom) of U."""
        flag = lambda keys: np.isin(self.keys, keys).astype(np.int8)
        return pd.DataFrame({'gene': np.asarray(self.genes, dtype=object)[self.keys >> POS_BITS],
                             'position': self.keys & ((1 << POS_BITS) - 1),
                             'nd': self.nd, 'observed_positive': flag(self.observed),
                             'observed_top': flag(self.top), 'observed_bottom': flag(self.bottom)})


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


def build_universe(genes, null_df, observed_df, background_file=None):
    if background_file is None:
        raise ValueError('Explicit tested positions are required')
    if len(genes) != len(set(genes)):
        raise ValueError('Duplicate background genes')
    for g in genes:
        identifier(g)
    gene_to_id = {g: i for i, g in enumerate(genes)}
    tested = read_tested(background_file, set(genes))
    keys = np.sort(np.concatenate([pair_keys(np.full(len(v), gene_to_id[g]), v)
                                   for g, v in tested.items()])) if tested else np.array([], dtype=np.int64)
    if observed_df.duplicated(['gene', 'position', 'side']).any():
        raise ValueError('Duplicate observed position side')
    inside, _ = _membership(observed_df, tested)
    if not inside.all():
        raise ValueError('Observed position outside declared tested universe / cleaned background')
    side_keys = {}
    for side in SIDES:
        rows = observed_df[observed_df.side == side]
        side_keys[side] = np.unique(pair_keys(rows.gene.map(gene_to_id), rows.position))
    if null_df.duplicated(['gene', 'position']).any():
        raise ValueError('Duplicate null positions')
    keep, _ = _membership(null_df, tested)
    rows = null_df[keep]
    null_keys = pair_keys(rows.gene.map(gene_to_id), rows.position)
    nd = np.zeros(len(keys), dtype=np.int64)
    nd[np.searchsorted(keys, null_keys)] = rows.nd.to_numpy()
    observed = np.union1d(*side_keys.values())
    metrics = {'n_genes': len(genes), 'n_universe': len(keys), 'tested_set': 'used',
               'n_tested_positions': len(keys), 'n_observed_positive': len(observed),
               'n_observed_outside_tested': 0, 'n_null_positive_all': len(null_df),
               'n_null_positive_in_universe': int(keep.sum()),
               'n_observed_not_scored_by_null': int((~np.isin(observed, null_keys)).sum())}
    return Universe(genes, gene_to_id, keys, nd, side_keys['top'], side_keys['bottom'], metrics)


def validate_null_inputs(observed, contract_path, null_path, roster_path, counts_path, blocks, output, background_positions):
    # Reuse the scientific/checksum checks already approved for SCORING.
    scoring = Path(__file__).resolve().parents[3]/'SCORING/local/src'
    sys.path.insert(0, str(scoring))
    from null_position_stats import validate_null, SCORE_RULE, read_tested as read_declared_tested
    import pyarrow.parquet as pq
    contract = json.loads(Path(contract_path).read_text())
    if contract.get('version') != 1 or contract.get('observed_sha256') != sha256_file(observed) or not contract.get('scientific'):
        raise ValueError('Observed provenance absent or checksum mismatch; regenerate CT_POSTPROC')
    roster = Path(roster_path).read_text().splitlines()
    counts = json.loads(Path(counts_path).read_text())
    if counts.get('score_rule') != SCORE_RULE or not blocks:
        raise ValueError('Canonical null blocks / score rule required')
    # the null must have been scored over the positional universe of this background.output
    validate_null([Path(b) for b in blocks], roster, pq.read_table(null_path), counts, output, contract,
                  read_declared_tested(background_positions))
    return {'N': len(roster), 'blocks': block_inventory(blocks), 'cycles': roster,
            'observed_contract': contract}


def position_calibration(observed, blocks, null_path, universe, N):
    """Independent consumer check of SCORING p/q on exactly U and the same cycles."""
    from null_position_stats import observed_best, count_block
    import pyarrow.parquet as pq
    genes, positions, scores = observed_best(Path(observed))
    classes = pq.read_table(null_path)
    ge = np.zeros(len(genes), dtype=np.int64)
    for block in blocks:
        counts = count_block(Path(block)/'positions.parquet', classes, np.array([]), genes, positions, scores)
        ge += counts['k_emp'].to_numpy()
    p = (ge + 1)/(N + 1)
    q = np.ones(len(p))
    if len(p):
        order = np.argsort(p, kind='stable')
        q[order] = np.minimum(1., np.minimum.accumulate((p[order]*len(universe.keys)/np.arange(1,len(p)+1))[::-1])[::-1])
    nd_lookup = dict(zip(universe.keys.tolist(), universe.nd.tolist()))
    keys = pair_keys([universe.gene_to_id[g] for g in genes], positions)
    nd = np.array([nd_lookup[int(key)] for key in keys],dtype=np.int64)
    if np.any(ge > nd) or np.any(nd > N):
        raise ValueError('Invalid positional calibration counts')
    return pd.DataFrame(dict(Gene=genes, Position=positions, CAAS_score=scores,
                             **{'p.emp':p,'p.adj_bh':q,'nd':nd,'k_emp':ge}))


def manifest(universe, inputs):
    """Manifest dict: SHA-256 of every input file, counts and the tested-set decision."""
    return {'component': 'universe',
            'contract_version': 2,
            'genes': universe.genes,
            'inputs': {k: {'path': os.path.abspath(v), 'sha256': sha256_file(v)} for k, v in inputs.items() if v},
            'gene_list_sha256': sha256_file(inputs['background_genes']),
            'gene_set_sha256': gene_list_digest(universe.genes),
            'metrics': universe.metrics}


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--background-genes', required=True, help='gene list, one gene per line')
    p.add_argument('--background-positions', required=True, help='caastools background.output (gene<TAB>tested positions)')
    p.add_argument('--null-positions', required=True, help='null_positions.parquet')
    p.add_argument('--observed', required=True, help='filtered_discovery.tsv or position_scores.tsv')
    p.add_argument('--blocks', nargs='+', required=True)
    p.add_argument('--roster', required=True)
    p.add_argument('--counts', required=True)
    p.add_argument('--observed-contract', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--log-level', default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'])
    a = p.parse_args()
    logging.basicConfig(level=getattr(logging, a.log_level), format='%(asctime)s - %(levelname)s - %(message)s')

    genes = load_universe_genes(a.background_genes)
    null_df = load_null_positive(a.null_positions)
    null = validate_null_inputs(a.observed, a.observed_contract, a.null_positions, a.roster, a.counts, a.blocks, a.output_dir, a.background_positions)
    if set(null['observed_contract'].get('excluded_genes', [])) & set(genes):
        raise ValueError('Cleaned background includes explicitly excluded observed genes')
    universe = build_universe(genes, null_df, load_observed_positive(a.observed), a.background_positions)

    os.makedirs(a.output_dir, exist_ok=True)
    universe.table().to_csv(os.path.join(a.output_dir, 'universe.tsv.gz'), sep='\t', index=False)
    inputs = {'background_genes': a.background_genes, 'background_positions': a.background_positions,
              'null_positions': a.null_positions, 'observed': a.observed}
    inputs.update(roster=a.roster, counts=a.counts, observed_contract=a.observed_contract)
    calibration = position_calibration(a.observed, a.blocks, a.null_positions, universe, null['N'])
    calibration.to_csv(Path(a.output_dir)/'position_calibration.tsv',sep='\t',index=False)
    meta = manifest(universe, inputs)
    meta['null'] = null
    meta['universe_sha256'] = sha256_file(Path(a.output_dir)/'universe.tsv.gz')
    finish_manifest(a.output_dir, 'universe_manifest.json', meta)
    logging.info(f"|U| = {universe.metrics['n_universe']} "
                 f"({universe.metrics['n_observed_positive']} observed positive positions)")


if __name__ == '__main__':
    main()
