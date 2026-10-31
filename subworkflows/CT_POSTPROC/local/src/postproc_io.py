"""Strict observed contract and quality-filter inputs for CT_POSTPROC."""
import csv
import hashlib
import json
import math
from pathlib import Path

import pandas as pd

SCORE_RULE = 'caas_pos_species_pss_v4'
TYPES = {'convergent', 'divergent', 'single', 'no_change', 'symplesiomorphic', 'cross_hypothesis'}


def read_table(path):
    with Path(path).open(newline='') as fh:
        first = fh.readline()
        sep = '\t' if first.count('\t') > first.count(',') else ','
        fh.seek(0)
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        if len(header) != len(set(header)):
            raise ValueError('Duplicate input columns')
        for row in reader:
            if len(row) != len(header):
                raise ValueError('Malformed input row')
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def positions(values):
    numbers = pd.to_numeric(values, errors='raise')
    if any(not math.isfinite(x) or x < 0 or x != int(x) for x in numbers):
        raise ValueError('Positions must be finite nonnegative integers')
    return numbers.astype('int64')


def coverage(df, gene='gene', pos='msa_pos'):
    return {str(g): {'positions': len(ps := sorted(set(map(int, values)))),
                    'sha256': hashlib.sha256(''.join(f'{p}\n' for p in ps).encode()).hexdigest()}
            for g, values in df.groupby(gene, sort=False)[pos]}


def load_observed(master, stamp):
    record = json.loads(Path(stamp).read_text())
    if (record.get('version') != 1 or record.get('complete') is not True
            or not isinstance(record.get('genes'), dict) or not isinstance(record.get('batches'), list)):
        raise ValueError('Invalid observed completion record')
    genes, settings = {}, None
    excluded = set(record.get('excluded_genes', []))
    keys = ('model', 'posterior_threshold', 'tree_sha256', 'taxid_sha256', 'weights_sha256', 'weighted')
    for batch in record['batches']:
        if (batch.get('version') != 1 or batch.get('complete') is not True
                or batch.get('score_rule') != SCORE_RULE or not isinstance(batch.get('genes'), dict)):
            raise ValueError('Incomplete or incompatible observed batch')
        threshold = batch.get('posterior_threshold')
        if (not isinstance(threshold, (int, float)) or isinstance(threshold, bool)
                or not math.isfinite(threshold) or not 0 <= threshold <= 1
                or not batch.get('model') or not batch.get('tree_sha256')
                or not isinstance(batch.get('weighted'), bool)):
            raise ValueError('Missing or invalid scientific settings')
        current = tuple(batch.get(key) for key in keys)
        if settings is not None and settings != current:
            raise ValueError('Inconsistent observed batch settings')
        settings = current
        if genes.keys() & batch['genes'].keys():
            raise ValueError('Overlapping observed batches')
        genes.update(batch['genes'])
        excluded.update(batch.get('excluded_genes', []))
    if genes != record['genes'] or genes.keys() & excluded:
        raise ValueError('Inconsistent observed gene coverage or exclusions')
    df = read_table(master)
    if {'Gene', 'Position'} & set(df):
        raise ValueError('Conflicting structural column aliases')
    if not {'gene', 'msa_pos', 'side', 'asr_path_score'} <= set(df):
        raise ValueError('Missing master keys or score')
    if any(not g or g != g.strip() for g in df.gene):
        raise ValueError('Invalid master gene')
    df['msa_pos'] = positions(df.msa_pos)
    scores = pd.to_numeric(df.asr_path_score, errors='raise')
    if any(not math.isfinite(x) or not 0 <= x <= 1 for x in scores):
        raise ValueError('Invalid asr_path_score')
    if not df.side.isin(['top', 'bottom', 'none']).all():
        raise ValueError('Invalid side')
    if ((df.side == 'none') & (scores != 0)).any() or ((df.side != 'none') & (scores <= 0)).any():
        raise ValueError('Score conflicts with side')
    if df.duplicated(['gene', 'msa_pos', 'side']).any():
        raise ValueError('Duplicate gene/position/side')
    none_keys = set(df.loc[df.side == 'none', ['gene', 'msa_pos']].itertuples(index=False, name=None))
    positive_keys = set(df.loc[df.side != 'none', ['gene', 'msa_pos']].itertuples(index=False, name=None))
    if none_keys & positive_keys:
        raise ValueError('Conflicting none-side rows')
    if 'agreement_ambiguous' in df and not df.agreement_ambiguous.str.lower().isin(['true','false','']).all():
        raise ValueError('Invalid agreement_ambiguous')
    if 'convergence_type' in df and not df.convergence_type.isin(TYPES | {''}).all():
        raise ValueError('Invalid convergence_type')
    for col, lower, upper, integer in [('derived_agreement', 0, 1, False), ('n_hypotheses', 1, math.inf, True)]:
        if col in df:
            vals = pd.to_numeric(df.loc[df[col] != '', col], errors='raise')
            if any(not math.isfinite(x) or not lower <= x <= upper or (integer and x != int(x)) for x in vals):
                raise ValueError(f'Invalid {col}')
    if coverage(df) != record['genes']:
        raise ValueError('Master position coverage differs from observed record')
    # Keep the original score text: serialization must not round scientific values.
    return df, record


def validate_cluster(discovery, clusters):
    if not {'Gene', 'Position', 'clustering_flag'} <= set(clusters):
        raise ValueError('Missing cluster columns')
    clusters = clusters.copy()
    clusters['Position'] = positions(clusters.Position)
    if not clusters.clustering_flag.isin(['Good', 'Discarded']).all():
        raise ValueError('Invalid cluster flag')
    if clusters.duplicated(['Gene', 'Position']).any():
        raise ValueError('Duplicate cluster positions')
    expected = set(discovery[['Gene', 'Position']].itertuples(index=False, name=None))
    actual = set(clusters[['Gene', 'Position']].itertuples(index=False, name=None))
    if expected != actual:
        raise ValueError('Cluster coverage differs from candidate positions')
    return clusters


def parameter_token(value):
    from decimal import Decimal
    text = format(Decimal(str(value)) * 100, 'f')
    return text.rstrip('0').rstrip('.') if '.' in text else text
