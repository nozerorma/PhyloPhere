#!/usr/bin/env python3
# posenrich_members.py — Memberships of the registry in their compact form, and their expansion.
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""The registry stores memberships by what they are, not as one row per (term, position):

  attribute collections  one value per position, a column of positions.tsv:
                           ucr_core, ucr_flank        0/1
                           fubar_call                 positive / purifying / neutral / '' (untested)
                           fade_top_sig, fade_bottom_sig  0/1
                           window                     genomic window name or ''
                         collections ucr_core, ucr_flank, fubar_positive, fubar_purifying, fubar_tested, fade_top_sig,
                         fade_bottom_sig and genomic_windows are derived from them (single term 'global', except the
                         windows, one term per window)
  region collections     regions.tsv: one record per run of consecutive CAAS columns of a term in a gene (Pfam domains
                         and clans); the positions of a run are start..end
  site collections       sites.tsv: one row per (term, gene, position) (custom lists)

`iter_members` yields one collection at a time and `read_members` concatenates them, so a consumer never holds the
memberships of every collection at once unless it asks for them.
"""

import numpy as np
import pandas as pd
from posenrich_contract import read_table

MEMBER_COLUMNS = ['collection_id', 'term_id', 'term_label', 'Gene', 'Position']
COVERAGE_COLLECTIONS = ('fubar_tested', 'genomic_windows')   # describe coverage, not annotation

# collection_id -> (column of positions.tsv, test of the column, term label)
ATTRIBUTE_COLLECTIONS = {
    'ucr_core': ('ucr_core', lambda c: c == '1', 'UCR_core'),
    'ucr_flank': ('ucr_flank', lambda c: c == '1', 'UCR_flank'),
    'fubar_positive': ('fubar_call', lambda c: c == 'positive', 'FUBAR_positive'),
    'fubar_purifying': ('fubar_call', lambda c: c == 'purifying', 'FUBAR_purifying'),
    'fubar_tested': ('fubar_call', lambda c: c != '', 'FUBAR_tested'),
    'fade_top_sig': ('fade_top_sig', lambda c: c == '1', 'fade_top_sig'),
    'fade_bottom_sig': ('fade_bottom_sig', lambda c: c == '1', 'fade_bottom_sig'),
    'genomic_windows': ('window', lambda c: c != '', None),
}
ATTRIBUTE_COLUMNS = ['ucr_core', 'custom_terms', 'ucr_flank', 'ucr_region', 'ucr_variability', 'fubar_call',
                     'fade_top_sig', 'fade_bottom_sig', 'window']
REGION_COLUMNS = ['collection_id', 'term_id', 'term_label', 'Gene', 'start', 'end',
                  'source_record_id', 'mapping_status', 'detail']
SITE_COLUMNS = ['collection_id', 'term_id', 'term_label', 'Gene', 'Position',
                'source_record_id', 'mapping_status', 'detail']


def attribute_members(positions, collection_id):
    """Members of an attribute collection from a positions table with its attribute columns (strings)."""
    column, test, label = ATTRIBUTE_COLLECTIONS[collection_id]
    values = positions[column].fillna('').astype(str)
    hit = values.map(test).to_numpy(dtype=bool)
    sub = positions.loc[hit, ['Gene', 'Position']]
    if collection_id == 'genomic_windows':
        term = values[hit].to_numpy()
        labels = term
    else:
        term = np.full(len(sub), 'global', dtype=object)
        labels = np.full(len(sub), label, dtype=object)
    out = pd.DataFrame({'collection_id': collection_id, 'term_id': term, 'term_label': labels,
                        'Gene': sub['Gene'].to_numpy(), 'Position': sub['Position'].astype(np.int64).to_numpy()})
    return out


def expand_runs(regions):
    """One row per (record, position) of the runs start..end, without repeats of (term, gene, position)."""
    if regions.empty:
        return pd.DataFrame(columns=MEMBER_COLUMNS)
    start = regions['start'].astype(np.int64).to_numpy()
    lengths = regions['end'].astype(np.int64).to_numpy() - start + 1
    which = np.repeat(np.arange(len(regions)), lengths)
    offset = np.arange(lengths.sum()) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    out = regions[['collection_id', 'term_id', 'term_label', 'Gene']].iloc[which].reset_index(drop=True)
    out['Position'] = start[which] + offset
    return out.drop_duplicates(['collection_id', 'term_id', 'Gene', 'Position']).reset_index(drop=True)


def iter_members(positions, regions, sites, collections=None):
    """Yield (collection_id, DataFrame[MEMBER_COLUMNS]) for each collection present, one at a time.

    positions: positions table with the attribute columns (None when no attribute collection is wanted),
    regions / sites: the region and site tables (or None). collections: ids wanted, None for all.
    """
    wanted = None if collections is None else set(collections)
    if positions is not None:
        for cid in ATTRIBUTE_COLLECTIONS:
            if (wanted is None or cid in wanted) and ATTRIBUTE_COLUMNS[0] in positions.columns:
                frame = attribute_members(positions, cid)
                if len(frame):
                    yield cid, frame
    if regions is not None and len(regions):
        for cid, sub in regions.groupby('collection_id', sort=True):
            if wanted is None or cid in wanted:
                yield cid, expand_runs(sub)
    if sites is not None and len(sites):
        for cid, sub in sites.groupby('collection_id', sort=True):
            if wanted is None or cid in wanted:
                frame = sub[MEMBER_COLUMNS].copy()
                frame['Position'] = frame['Position'].astype(np.int64)
                yield cid, frame.drop_duplicates(['collection_id', 'term_id', 'Gene', 'Position']).reset_index(drop=True)


def read_tables(registry_dir, collections=None):
    """(positions, regions, sites) of a registry directory; only what the wanted collections need is read."""
    import os
    wanted = None if collections is None else set(collections)
    need_attr = wanted is None or bool(wanted & set(ATTRIBUTE_COLLECTIONS))
    positions = None
    if need_attr:
        usecols = ['Gene', 'Position'] + ATTRIBUTE_COLUMNS
        positions = read_table(os.path.join(registry_dir, 'positions.tsv'))
        positions = positions[[c for c in usecols if c in positions]]
    regions = read_table(os.path.join(registry_dir, 'regions.tsv'))
    sites = read_table(os.path.join(registry_dir, 'sites.tsv'))
    return positions, regions, sites


def read_members(registry_dir, collections=None, genes=None):
    """All members of the wanted collections as one DataFrame[MEMBER_COLUMNS], optionally of some genes only."""
    positions, regions, sites = read_tables(registry_dir, collections)
    if genes is not None:
        genes = set(genes)
        if positions is not None:
            positions = positions[positions['Gene'].isin(genes)]
        if regions is not None:
            regions = regions[regions['Gene'].isin(genes)]
        if sites is not None:
            sites = sites[sites['Gene'].isin(genes)]
    parts = [frame for _, frame in iter_members(positions, regions, sites, collections)]
    if not parts:
        return pd.DataFrame(columns=MEMBER_COLUMNS)
    return pd.concat(parts, ignore_index=True)


def member_records(registry_dir, collections=None, genes=None):
    """Expanded members as dict rows with string values: collection_id, term_id, term_label, Gene, Position."""
    frame = read_members(registry_dir, collections, genes)
    for cid, term, label, gene, pos in zip(frame['collection_id'], frame['term_id'], frame['term_label'],
                                           frame['Gene'], frame['Position']):
        yield {'collection_id': cid, 'term_id': term, 'term_label': label, 'Gene': gene, 'Position': str(int(pos))}
