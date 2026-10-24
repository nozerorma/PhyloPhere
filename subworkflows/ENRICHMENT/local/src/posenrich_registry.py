#!/usr/bin/env python3
# posenrich_registry.py — Common registry of positions, collections and coverage for POSENRICH modular.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Build the normalized position registry shared by every POSENRICH component.

WP10 P1 (docs/fases/WP10_posenrich_modular.md). The registry keeps normalized
tables; GMT is only a compatible export (--gmt_dir), never the sole record.

Inputs (all optional except --output_dir; NO_* sentinels and missing paths mean absent):
  --position_scores        SCORING's position_scores.tsv (one row per Gene, Position, side
                           with CAAS_score > 0; p.emp shared between sides, p.adj_bh, nd, ...)
  --postproc               filtered_discovery.tsv: positions evaluated by CT (a row may
                           have no positive score: evaluated without detection)
  --null_universe          null_universe_summary.tsv: N cycles, |U|, observed/unseen counts
  --map_dir                directory with <GENE>*.map.tsv (CAAS trimmed alignment)
  --cleaned_background     universe of genes (one per line)
  --gene_ensembl_file      gene/human_protein_id table (trazability only)
  --domain_variability_file, --ucr_positions_file, --fubar_sites_file,
  --fade_sites_top_file, --fade_sites_bottom_file, --custom_marker_file
                           the same collection sources as build_position_gmt.py

Outputs (in --output_dir):
  positions.tsv            one row per (Gene, Position) of the universe: evaluation status,
                           coordinates (raw codon column, reference residue, genomic) and MAP
  position_evidence.tsv    one row per (Gene, Position, side): CAAS_score, p.emp/p.adj_bh
                           (shared p, never two directional tests), nd, flags, N cycles
  collections.tsv          one row per (collection_id, term_id): type, coordinate space,
                           source hash, sizes, compatible analyses, within-protein contrast
  members.tsv              one row per (collection_id, term_id, Gene, Position): several
                           evidences may support one member without multiplying the row
  coverage.tsv             per gene/collection/universe metrics with explicit denominators
  pending_mappings.tsv     original records that could not be mapped; never turned into
                           negatives
  manifest.json            source inventory with SHA-256, included/empty/oversized/
                           unmapped/incompatible collections and audited losses
"""

import argparse
import csv
import hashlib
import json
import os
import sys

import pandas as pd

import build_position_gmt as bpg

CONTRACT_VERSION = 1

# collection_id -> (term_type, coordinate_space, analyses, informative_within_protein, source key)
COLLECTION_META = {
    'pfam_domains':       ('region', 'CAAS_0based', 'localization;global;proximity_protein', 1, 'domain_variability_file'),
    'pfam_clans':         ('region', 'CAAS_0based', 'localization;global;proximity_protein', 1, 'domain_variability_file'),
    'ucr_core':           ('region', 'CAAS_0based', 'localization;global;proximity_protein', 1, 'ucr_positions_file'),
    'ucr_flank':          ('region', 'CAAS_0based', 'localization;global;proximity_protein', 1, 'ucr_positions_file'),
    'fubar_positive':     ('site',   'CAAS_0based', 'localization;global;proximity_protein', 1, 'fubar_sites_file'),
    'fubar_purifying':    ('site',   'CAAS_0based', 'localization;global;proximity_protein', 1, 'fubar_sites_file'),
    'fubar_tested':       ('coverage', 'CAAS_0based', 'characterization', 0, 'fubar_sites_file'),
    'fade_top_sig':       ('site',   'CAAS_0based', 'localization;global;proximity_protein', 1, 'fade_sites_top_file'),
    'fade_bottom_sig':    ('site',   'CAAS_0based', 'localization;global;proximity_protein', 1, 'fade_sites_bottom_file'),
    'genomic_windows_1mb': ('window', 'hg38_1mb_bins', 'proximity_genomic', 0, 'map_dir'),
    'custom':             ('site',   'CAAS_0based', 'localization;global;proximity_protein', 1, 'custom_marker_file'),
}

# Global characterization layers (single-term collections) -> GMT layer name and description.
LAYER_NAMES = {
    'ucr_core': 'UCR_core', 'ucr_flank': 'UCR_flank',
    'fubar_positive': 'FUBAR_positive', 'fubar_purifying': 'FUBAR_purifying',
    'fade_top_sig': 'FADE_top_sig', 'fade_bottom_sig': 'FADE_bottom_sig',
}

# Coverage collections without a GMT layer: fixed description.
COVERAGE_LAYERS = {'fubar_tested': 'FUBAR tested positions (coverage)'}


def member_field(m, key):
    return m.get(key, '')

SOURCE_KEYS = ['position_scores', 'postproc', 'null_universe', 'map_dir', 'cleaned_background',
               'gene_ensembl_file', 'domain_variability_file', 'ucr_positions_file',
               'fubar_sites_file', 'fade_sites_top_file', 'fade_sites_bottom_file',
               'custom_marker_file']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for key in SOURCE_KEYS:
        p.add_argument(f'--{key}', required=False, default=None)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--gmt_dir', required=False, default=None,
                   help='Optional directory for the compatible GMT export of the registry')
    return p.parse_args()


def clean_optional(value):
    """NO_* sentinel, 'None' or empty means absent; resolve .gz siblings."""
    if value is None or value == '' or os.path.basename(str(value)).startswith('NO_'):
        return None
    return bpg.resolve_path(str(value))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


def load_position_keys(path, active_genes):
    """(gene, position) keys of position_scores.tsv plus the audited losses.

    Returns (keys, n_rows, n_outside_universe). Only CAAS_score > 0 rows exist in
    this table by contract.
    """
    keys, n_rows, n_out = set(), 0, 0
    if not path:
        return keys, n_rows, n_out
    for chunk in pd.read_csv(path, sep='\t', usecols=['Gene', 'Position'], dtype=str,
                             chunksize=500_000):
        for gene, pos in chunk.itertuples(index=False):
            n_rows += 1
            if gene not in active_genes:
                n_out += 1
                continue
            keys.add((gene, int(pos)))
    return keys, n_rows, n_out


def load_discovery_keys(path, active_genes):
    """(gene, position) keys of filtered_discovery.tsv (evaluated by CT)."""
    keys, n_rows, n_out = set(), 0, 0
    if not path:
        return keys, n_rows, n_out
    for chunk in pd.read_csv(path, sep='\t', usecols=['Gene', 'Position'], dtype=str,
                             chunksize=500_000):
        for gene, pos in chunk.itertuples(index=False):
            n_rows += 1
            if gene not in active_genes:
                n_out += 1
                continue
            keys.add((gene, int(pos)))
    return keys, n_rows, n_out


def load_null_universe(path):
    """{metric: value} of null_universe_summary.tsv (rows with no cycle count)."""
    if not path:
        return {}
    df = read_tsv(path)
    out = {}
    for row in df.itertuples(index=False):
        if str(row.cycles) in ('NA', ''):
            try:
                out[row.metric] = float(row.value)
            except (ValueError, TypeError):
                out[row.metric] = row.value
    return out


def build_positions(active_genes, map_cache, discovery_keys, score_keys):
    """positions.tsv rows: the union of MAP columns, CT detections and positive scores."""
    rows = {}

    def entry(gene, pos):
        return rows.setdefault((gene, pos), {
            'Gene': gene, 'Position': pos, 'in_map': 0, 'in_discovery': 0,
            'in_scores_positive': 0, 'ori_codon_col': 'NA', 'hg38_aa_pos': 'NA',
            'hg38_nt_coord': 'NA', 'strand': 'NA', 'map_file': 'NA'})

    for gene in sorted(active_genes):
        cache = map_cache.get(gene)
        if not cache:
            continue
        raw_to_col = cache['raw_to_col']
        col_to_raw = {c: r for r, c in raw_to_col.items()}
        col_to_residue = {c: r for r, c in cache['residue_to_col'].items()}
        for col in cache['selected_cols']:
            e = entry(gene, col)
            e['in_map'] = 1
            e['ori_codon_col'] = str(col_to_raw.get(col, 'NA'))
            e['hg38_aa_pos'] = str(col_to_residue.get(col, 'NA'))
            e['hg38_nt_coord'] = cache['col_to_genomic'].get(col, 'NA')
            e['strand'] = cache['strand']
            e['map_file'] = os.path.basename(cache.get('map_file', 'NA'))

    for gene, pos in discovery_keys:
        entry(gene, pos)['in_discovery'] = 1
    for gene, pos in score_keys:
        entry(gene, pos)['in_scores_positive'] = 1

    out = []
    for (gene, pos), e in sorted(rows.items()):
        if e['in_scores_positive']:
            status = 'detected_positive'
        elif e['in_discovery']:
            status = 'evaluated_no_detection'
        elif e['in_map']:
            status = 'map_only_not_in_discovery'
        else:
            status = 'not_evaluable'
        cause = '' if e['in_map'] else 'no_map'
        out.append({'Gene': gene, 'Position': pos, 'evaluation_status': status,
                    'exclusion_cause': cause, 'in_map': e['in_map'],
                    'in_discovery': e['in_discovery'],
                    'in_scores_positive': e['in_scores_positive'],
                    'ori_codon_col': e['ori_codon_col'], 'hg38_aa_pos': e['hg38_aa_pos'],
                    'hg38_nt_coord': e['hg38_nt_coord'], 'strand': e['strand'],
                    'map_file': e['map_file']})
    return pd.DataFrame(out, columns=['Gene', 'Position', 'evaluation_status', 'exclusion_cause',
                                      'in_map', 'in_discovery', 'in_scores_positive',
                                      'ori_codon_col', 'hg38_aa_pos', 'hg38_nt_coord',
                                      'strand', 'map_file'])


def build_evidence(position_scores_path, active_genes, n_cycles, postproc_path=None):
    """Canonical positive scores plus audited discovery rows without a score.

    core.scores.position_score defines nonpositive values as unavailable, not zero.
    Keep the raw path score separately; never impute a canonical score or a p-value.
    """
    cols = ['Gene', 'Position', 'side', 'CAAS_score', 'p.emp', 'p.adj_bh', 'p.emp_fact',
            'p.adj_bh_fact', 'nd', 'null_unseen', 'caas_nominal', 'caas_adjusted']
    extra = ['p_shared', 'N_cycles', 'evidence_status', 'score_status', 'raw_asr_path_score']
    df = read_tsv(position_scores_path) if position_scores_path else pd.DataFrame(columns=cols)
    df = df.reindex(columns=cols).fillna('NA')
    df = df[df['Gene'].isin(active_genes)].copy()
    df['Position'] = df['Position'].astype(int)
    if df.duplicated(['Gene', 'Position', 'side']).any():
        raise ValueError('duplicate position score side; one row per Gene/Position/side required')
    scores = pd.to_numeric(df['CAAS_score'], errors='coerce')
    if ((scores <= 0) | ~scores.notna()).any():
        raise ValueError('position_scores must contain canonical positive CAAS scores')
    for key, sub in df.groupby(['Gene', 'Position']):
        for col in ('p.emp', 'p.adj_bh'):
            if sub[col].nunique() > 1:
                raise ValueError(f'inconsistent shared {col} for {key}')
    df['score_status'] = 'canonical_positive'
    df['raw_asr_path_score'] = 'NA'
    rows = df.to_dict('records')
    by_key = {(row['Gene'], int(row['Position']), row['side']): row for row in rows}
    if postproc_path:
        discovery = read_tsv(postproc_path)
        discovery = discovery[discovery['Gene'].isin(active_genes)]
        if discovery.duplicated(['Gene', 'Position', 'side']).any():
            raise ValueError('duplicate discovery side; cannot resolve canonical score')
        for original in discovery.to_dict('records'):
            key = (original['Gene'], int(original['Position']), original['side'])
            raw = original.get('asr_path_score', 'NA')
            number = pd.to_numeric(raw, errors='coerce')
            if key in by_key:
                row = by_key[key]
                if pd.notna(number) and number > 0 and abs(float(row['CAAS_score']) - number) > 1e-12:
                    raise ValueError(f'discovery and canonical score disagree for {key}')
                row['raw_asr_path_score'] = raw
                continue
            row = {c: 'NA' for c in cols}
            row.update(Gene=key[0], Position=key[1], side=key[2], raw_asr_path_score=raw)
            if original['side'] not in ('top', 'bottom') or pd.isna(number):
                row['score_status'] = 'not_evaluable'
            elif number <= 0:
                row['score_status'] = 'evaluated_no_canonical_score'
            else:
                # Same identity rule as core.scores.position_score, only for positive sides.
                row['CAAS_score'] = repr(float(number))
                row['score_status'] = 'canonical_positive_from_discovery'
            rows.append(row)
    for row in rows:
        row['p_shared'] = 'TRUE'
        row['N_cycles'] = 'NA' if n_cycles is None else str(int(n_cycles))
        row['evidence_status'] = 'evaluated' if row['p.emp'] not in ('NA', '') else 'p_unavailable'
    # A shared p is copied to any additional side of the same biological position.
    shared = {(row['Gene'], int(row['Position'])): {c: row[c] for c in cols[4:]}
              for row in rows if row['p.emp'] not in ('NA', '')}
    for row in rows:
        values = shared.get((row['Gene'], int(row['Position'])))
        if values:
            row.update(values)
            row['evidence_status'] = 'evaluated'
    return pd.DataFrame(rows, columns=cols + extra).sort_values(['Gene', 'Position', 'side'])


def collect_all_members(args, map_cache, active_genes, pending):
    """Run every collection builder, gathering members, terms and char layers."""
    result = {}
    result['pfam_domains'], result['pfam_clans'], pfam_char, m = bpg.collect_pfam_members(
        args.domain_variability_file, map_cache, active_genes, pending=pending)
    result['_pfam_char'] = pfam_char
    result['members'] = list(m)
    result['ucr_core_cols'], result['ucr_flank_cols'], ucr_region_char, ucr_variability_char, m = \
        bpg.collect_ucr_members(args.ucr_positions_file, map_cache, active_genes, pending=pending)
    result['_ucr_region_char'] = ucr_region_char
    result['_ucr_variability_char'] = ucr_variability_char
    result['members'] += m
    result['fubar_positive_cols'], result['fubar_purifying_cols'], fubar_char, m = \
        bpg.collect_fubar_members(args.fubar_sites_file, map_cache, active_genes, pending=pending)
    result['_fubar_char'] = fubar_char
    result['members'] += m
    result['fade_top_cols'], m = bpg.collect_fade_sites(
        args.fade_sites_top_file, map_cache, active_genes, 'fade_top_sig', pending=pending)
    result['members'] += m
    result['fade_bottom_cols'], m = bpg.collect_fade_sites(
        args.fade_sites_bottom_file, map_cache, active_genes, 'fade_bottom_sig', pending=pending)
    result['members'] += m
    result['genomic_windows_1mb'], m = bpg.collect_genomic_windows(map_cache)
    result['members'] += m
    custom_error = None
    try:
        result['custom'], m = bpg.collect_custom_members(
            args.custom_marker_file, map_cache, active_genes, strict=False, pending=pending)
        result['members'] += m
    except ValueError as e:
        custom_error = str(e)
    result['_custom_error'] = custom_error
    result['_char_layers'] = bpg.build_char_layers(
        result['ucr_core_cols'], result['ucr_flank_cols'],
        result['fubar_positive_cols'], result['fubar_purifying_cols'],
        result['fade_top_cols'], result['fade_bottom_cols'])
    return result


def build_members_df(members):
    """members.tsv: one row per (collection_id, term_id, Gene, Position).

    Several evidences may support one member without multiplying the row: their
    source_record_id values are joined with ';'.
    """
    agg = {}
    for m in members:
        key = (m['collection_id'], m['term_id'], m['gene'], m['position'])
        cur = agg.get(key)
        if cur is None:
            agg[key] = {'records': [m['source_record_id']],
                        'labels': [member_field(m, 'term_label')],
                        'details': [member_field(m, 'detail')],
                        'status': m['mapping_status']}
        else:
            if m['source_record_id'] not in cur['records']:
                cur['records'].append(m['source_record_id'])
            if member_field(m, 'term_label') not in cur['labels']:
                cur['labels'].append(member_field(m, 'term_label'))
            if member_field(m, 'detail') and member_field(m, 'detail') not in cur['details']:
                cur['details'].append(member_field(m, 'detail'))
            if m['mapping_status'] == 'direct_caas':
                cur['status'] = 'direct_caas'
    rows = [{'collection_id': c, 'term_id': t, 'term_label': ';'.join(sorted(set(v['labels']))),
             'Gene': g, 'Position': p,
             'source_record_id': ';'.join(sorted(v['records'])),
             'mapping_status': v['status'],
             'detail': ';'.join(sorted(set(v['details'])))}
            for (c, t, g, p), v in sorted(agg.items())]
    return pd.DataFrame(rows, columns=['collection_id', 'term_id', 'term_label', 'Gene', 'Position',
                                       'source_record_id', 'mapping_status', 'detail'])


def build_collections_df(coll, members_df, sources):
    """collections.tsv: one row per (collection_id, term_id)."""
    rows = []
    grouped = members_df.groupby(['collection_id', 'term_id'])
    # Terms present in the GMT-ready dicts but without members (omitted from GMT).
    term_desc = {}
    for cid in ('pfam_domains', 'pfam_clans', 'genomic_windows_1mb', 'custom'):
        for term, (desc, _members) in coll.get(cid, {}).items():
            term_desc[(cid, term)] = desc
    for cid, layer in LAYER_NAMES.items():
        source_key = COLLECTION_META[cid][4]
        if sources.get(source_key, 'NA') == 'NA':
            continue  # absent source: the collection is missing, not an empty term
        desc = coll['_char_layers'].get(layer, ('', []))[0]
        term_desc[(cid, 'global')] = desc
    for cid, desc in COVERAGE_LAYERS.items():
        source_key = COLLECTION_META[cid][4]
        if sources.get(source_key, 'NA') == 'NA':
            continue
        term_desc[(cid, 'global')] = desc
    if coll.get('_custom_error') and not any(cid == 'custom' for cid, _ in term_desc):
        # The custom source was rejected: report the collection, not a silent absence.
        term_desc[('custom', 'NA')] = 'Custom annotation term'
    for (cid, term), desc in sorted(term_desc.items()):
        meta = COLLECTION_META[cid]
        sub = grouped.get_group((cid, term)) if (cid, term) in grouped.groups else None
        n_members = 0 if sub is None else len(sub)
        n_genes = 0 if sub is None else sub['Gene'].nunique()
        status = 'included' if n_members else 'empty'
        if cid == 'custom' and coll.get('_custom_error'):
            status = 'incompatible'
        rows.append({'collection_id': cid, 'term_id': term, 'description': desc,
                     'term_type': meta[0], 'coordinate_space': meta[1],
                     'source_file': sources.get(meta[4], 'NA'),
                     'source_sha256': sources.get(meta[4] + '_sha256', 'NA'),
                     'n_members': n_members, 'n_genes': n_genes,
                     'analyses': meta[2], 'informative_within_protein': meta[3],
                     'status': status})
    return pd.DataFrame(rows, columns=['collection_id', 'term_id', 'description', 'term_type',
                                       'coordinate_space', 'source_file', 'source_sha256',
                                       'n_members', 'n_genes', 'analyses',
                                       'informative_within_protein', 'status'])


def build_coverage_df(positions, members_df, collections_df, null_universe, universe_genes):
    """coverage.tsv: metrics with explicit denominators and their limitations."""
    rows = []

    def fmt(v):
        if isinstance(v, float) and v.is_integer():
            return str(int(v))
        return str(v)

    def add(scope, ident, metric, value, denominator, denom_def, source, limitation):
        rows.append({'scope': scope, 'id': ident, 'metric': metric, 'value': fmt(value),
                     'denominator': fmt(denominator), 'denominator_definition': denom_def,
                     'source': source, 'limitation': limitation})

    n_map = int(positions['in_map'].sum()) if not positions.empty else 0
    n_disc = int(positions['in_discovery'].sum()) if not positions.empty else 0
    n_pos = int(positions['in_scores_positive'].sum()) if not positions.empty else 0
    denom_def = 'suma de columnas seleccionadas del MAP sobre los genes del universo (posiciones anotables)'
    add('universe', '*', 'n_genes', len(universe_genes), 'NA', 'genes del universo (cleaned_background)', 'cleaned_background', '')
    add('universe', '*', 'n_positions_total', len(positions), 'NA', 'posiciones del registro', 'registry', '')
    add('universe', '*', 'n_positions_annotable', n_map, 'NA', denom_def, 'map_dir',
        'desconocida sin MAP; no igualar a cobertura experimental')
    add('universe', '*', 'n_positions_evaluated_discovery', n_disc, 'NA', 'posiciones detectadas por CT', 'filtered_discovery.tsv',
        'evaluada por CT; no es el universo completo de posiciones')
    add('universe', '*', 'n_positions_positive', n_pos, 'NA', 'posiciones con CAAS_score > 0', 'position_scores.tsv',
        'position_scores.tsv no contiene todas las posiciones evaluadas')
    annotated = set()
    if not members_df.empty:
        annotated = set(zip(members_df['Gene'], members_df['Position']))
    n_annot = sum(1 for g, p in zip(positions['Gene'], positions['Position']) if (g, p) in annotated) \
        if not positions.empty else 0
    add('universe', '*', 'n_positions_annotated', n_annot, 'NA', 'posiciones en al menos una coleccion', 'members.tsv',
        'ausencia de anotacion no demuestra ausencia biologica')
    for metric in ('N_cycles', 'null_positive_positions_in_universe', 'observed_positive_positions',
                   'observed_unseen_positions', 'universe_size'):
        if metric in null_universe:
            add('universe', '*', metric, int(null_universe[metric]), 'NA',
                'metricas del nulo posicional canónico', 'null_universe_summary.tsv',
                'el universo del null no equivale al de coordenadas; no reconstruye scores por ciclo')
    if not positions.empty:
        for gene, sub in positions.groupby('Gene'):
            add('gene', gene, 'n_map_positions', int(sub['in_map'].sum()), 'NA', denom_def, 'map_dir', '')
            add('gene', gene, 'n_discovery_positions', int(sub['in_discovery'].sum()), 'NA',
                'posiciones detectadas por CT', 'filtered_discovery.tsv', '')
            add('gene', gene, 'n_positive_positions', int(sub['in_scores_positive'].sum()), 'NA',
                'posiciones con CAAS_score > 0', 'position_scores.tsv', '')
            n_ann = sum(1 for g, p in zip(sub['Gene'], sub['Position']) if (g, p) in annotated)
            add('gene', gene, 'n_annotated_positions', n_ann, 'NA', 'posiciones en al menos una coleccion', 'members.tsv', '')
    if not members_df.empty:
        for cid, sub in members_df.groupby('collection_id'):
            n_positions = len(set(zip(sub['Gene'], sub['Position'])))
            prop = (n_positions / n_map) if n_map else 'NA'
            add('collection', cid, 'n_terms', int(sub['term_id'].nunique()), 'NA', 'terminos con al menos un miembro', 'members.tsv', '')
            add('collection', cid, 'n_positions', n_positions, n_map, denom_def, 'members.tsv', '')
            add('collection', cid, 'coverage_proportion', prop, n_map, denom_def, 'members.tsv',
                'proporcion de posiciones anotables cubiertas por la coleccion')
            add('collection', cid, 'n_genes', int(sub['Gene'].nunique()), 'NA', 'genes con al menos un miembro', 'members.tsv', '')
    return pd.DataFrame(rows, columns=['scope', 'id', 'metric', 'value', 'denominator',
                                       'denominator_definition', 'source', 'limitation'])


def export_gmt(gmt_dir, coll):
    """Compatible GMT export from the registry terms (same files as build_position_gmt.py)."""
    os.makedirs(gmt_dir, exist_ok=True)
    bpg.write_gmt(os.path.join(gmt_dir, 'pfam_domains.gmt'), coll['pfam_domains'])
    bpg.write_gmt(os.path.join(gmt_dir, 'pfam_clans.gmt'), coll['pfam_clans'])
    bpg.write_gmt(os.path.join(gmt_dir, 'genomic_locations.gmt'), coll['genomic_windows_1mb'])
    if coll.get('custom'):
        bpg.write_gmt(os.path.join(gmt_dir, 'custom_features.gmt'), coll['custom'])
    bpg.write_gmt(os.path.join(gmt_dir, 'characterization_layers.tsv'), coll['_char_layers'])


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    for key in SOURCE_KEYS:
        setattr(args, key, clean_optional(getattr(args, key)))
    if args.map_dir and not os.path.isdir(args.map_dir):
        print(f"WARNING: --map_dir is not a directory (skipping coordinate mapping): {args.map_dir}",
              file=sys.stderr)
        args.map_dir = None

    sources = {}
    for key in SOURCE_KEYS:
        path = getattr(args, key)
        if path is None:
            sources[key] = 'NA'
            sources[key + '_sha256'] = 'NA'
            continue
        sources[key] = path
        sources[key + '_sha256'] = sha256_file(path) if os.path.isfile(path) else 'NA(directory)'

    universe_genes = bpg.load_universe_genes(args.cleaned_background)
    map_cache = bpg.build_map_cache(args.map_dir, universe_genes)
    active_genes = set(universe_genes) if universe_genes is not None else set(map_cache.keys())

    discovery_keys, n_disc_rows, n_disc_out = load_discovery_keys(args.postproc, active_genes)
    score_keys, n_score_rows, n_score_out = load_position_keys(args.position_scores, active_genes)
    null_universe = load_null_universe(args.null_universe)
    n_cycles = null_universe.get('N_cycles')

    pending = []
    coll = collect_all_members(args, map_cache, active_genes, pending)

    positions = build_positions(active_genes, map_cache, discovery_keys, score_keys)
    evidence = build_evidence(args.position_scores, active_genes, n_cycles, args.postproc)
    if not evidence.empty:
        positive = evidence[evidence['CAAS_score'] != 'NA']
        positive_keys = set(zip(positive['Gene'], positive['Position'].astype(int)))
        for idx, row in positions.iterrows():
            if (row['Gene'], row['Position']) in positive_keys:
                positions.loc[idx, 'evaluation_status'] = 'detected_positive'
        evaluable = evidence[evidence['score_status'] != 'not_evaluable']
        known = set(zip(evaluable['Gene'], evaluable['Position'].astype(int)))
        for idx, row in positions.iterrows():
            if row['in_discovery'] and not row['in_scores_positive'] and (row['Gene'], row['Position']) not in known:
                positions.loc[idx, 'evaluation_status'] = 'not_evaluable'
    members_df = build_members_df(coll['members'])
    collections_df = build_collections_df(coll, members_df, sources)
    coverage_df = build_coverage_df(positions, members_df, collections_df, null_universe,
                                    active_genes)
    pending_df = pd.DataFrame(pending, columns=['source', 'record_id', 'gene', 'raw_position', 'reason']) \
        if pending else pd.DataFrame(columns=['source', 'record_id', 'gene', 'raw_position', 'reason'])
    pending_df['candidates'] = 'NA'
    pending_df['note'] = 'no mapeado a coordenadas CAAS; no contar como negativo'
    pending_df = pending_df[['source', 'record_id', 'gene', 'raw_position', 'reason',
                             'candidates', 'note']]

    write_tsv(positions, os.path.join(args.output_dir, 'positions.tsv'))
    write_tsv(evidence, os.path.join(args.output_dir, 'position_evidence.tsv'))
    write_tsv(collections_df, os.path.join(args.output_dir, 'collections.tsv'))
    write_tsv(members_df, os.path.join(args.output_dir, 'members.tsv'))
    write_tsv(coverage_df, os.path.join(args.output_dir, 'coverage.tsv'))
    write_tsv(pending_df, os.path.join(args.output_dir, 'pending_mappings.tsv'))
    if args.gmt_dir:
        export_gmt(args.gmt_dir, coll)

    # Manifest: included/empty/oversized/unmapped/incompatible, with audited losses.
    status_counts = positions['evaluation_status'].value_counts().to_dict() if not positions.empty else {}
    collections_manifest = []
    for cid, meta in COLLECTION_META.items():
        sub = collections_df[collections_df['collection_id'] == cid] if not collections_df.empty else \
            pd.DataFrame(columns=['term_id', 'n_members', 'status'])
        source_present = sources.get(meta[4], 'NA') != 'NA'
        entry = {'collection_id': cid,
                 'status': ('missing' if not source_present else
                            'incompatible' if cid == 'custom' and coll.get('_custom_error') else
                            'included' if not sub.empty and (sub['n_members'] > 0).any() else 'empty'),
                 'n_terms': int(len(sub)), 'n_members': int(sub['n_members'].sum()) if len(sub) else 0,
                 'n_genes': int(members_df[members_df['collection_id'] == cid]['Gene'].nunique())
                 if not members_df.empty else 0,
                 'source': sources.get(meta[4], 'NA')}
        if cid == 'custom' and coll.get('_custom_error'):
            entry['note'] = coll['_custom_error']
        collections_manifest.append(entry)
    manifest = {
        'contract_version': CONTRACT_VERSION,
        'coordinate_space': 'CAAS_0based',
        'universe': {'source': sources.get('cleaned_background', 'NA'), 'n_genes': len(active_genes)},
        'sources': {key: {'path': sources.get(key, 'NA'), 'sha256': sources.get(key + '_sha256', 'NA')}
                    for key in SOURCE_KEYS},
        'positions': {'total': int(len(positions)), 'by_status': {k: int(v) for k, v in status_counts.items()}},
        'evidence': {'rows': int(len(evidence)),
                     'ranking_population': 'canonical_positive_detections_only',
                     'complete_evaluable_universe': False,
                     'nonpositive_rule': 'canonical score unavailable; never imputed as zero',
                     'null_cycles': 'NA' if n_cycles is None else int(n_cycles),
                     'p_available': int((evidence['evidence_status'] == 'evaluated').sum())
                     if not evidence.empty else 0},
        'collections': collections_manifest,
        'pending_mappings': int(len(pending_df)),
        'audited_losses': {
            'discovery_rows_outside_universe': n_disc_out,
            'score_rows_outside_universe': n_score_out,
            'discovery_rows': n_disc_rows,
            'score_rows': n_score_rows,
        },
        'gmt_export': os.path.abspath(args.gmt_dir) if args.gmt_dir else None,
    }
    with open(os.path.join(args.output_dir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write('\n')

    print(f"registry: {len(positions)} positions ({len(active_genes)} genes), "
          f"{len(members_df)} member rows, {len(collections_df)} terms, "
          f"{len(pending_df)} pending mappings -> {args.output_dir}")


if __name__ == '__main__':
    main()
