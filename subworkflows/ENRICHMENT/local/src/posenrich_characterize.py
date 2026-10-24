#!/usr/bin/env python3
# posenrich_characterize.py — General characterization of prioritized positions (POSENRICH component 1).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Characterize the prioritized positions and the collections, without inference.

WP10 P2 (docs/fases/WP10_posenrich_modular.md). Consumes the common registry
(posenrich_registry.py) and produces the characterization tables. This component
does not run any enrichment test and assigns no enrichment significance to a
mere annotation. It works with or without a null, without hits and without
annotations: every table reports explicit states (executed without evidence,
not available, not executed).

Inputs:
  --registry_dir   directory with positions.tsv, position_evidence.tsv, members.tsv,
                   collections.tsv, coverage.tsv, manifest.json
  --output_dir     destination directory (required)

Outputs (in --output_dir):
  position_characterization.tsv  one row per annotated or scored (Gene, Position):
                                 annotation layers + ranking columns. Single owner of
                                 this table for reports 14/15.
  collection_summary.tsv         per collection/term: sizes, exact crosses with the
                                 ranking (counts and proportions, denominators shown)
  ranking_summary.tsv            score distribution and nominal/adjusted counts per
                                 ranking (global deduplicated, top, bottom); no
                                 selection percentiles
  universe_summary.tsv           evaluation states of the position universe
  characterization_manifest.json what was produced, from which inputs, and which
                                 components are not executed
"""

import argparse
import csv
import hashlib
import json
import os

import pandas as pd

REGISTRY_FILES = ['positions.tsv', 'position_evidence.tsv', 'members.tsv',
                  'collections.tsv', 'coverage.tsv', 'manifest.json']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--registry_dir', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--nominal_threshold', type=float, default=0.05,
                   help='Nominal p.emp threshold, only when flags are absent from the evidence')
    p.add_argument('--fdr_threshold', type=float, default=0.05,
                   help='BH threshold, only when flags are absent from the evidence')
    return p.parse_args()


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_detail(detail):
    """'k1=v1;k2=v2' -> dict."""
    out = {}
    for part in str(detail).split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            out[k] = v
    return out


def position_scores_by_key(evidence):
    """(gene, pos) -> {'score': best-side CAAS_score, 'sides': 'top;bottom', p.emp, ...}.

    CAAS_score stays directional in the evidence rows; the characterization table
    carries the best-side score plus the side list, while p.emp/p.adj_bh are the
    shared max-over-sides values (one test per position, not two).
    """
    out = {}
    for row in evidence.to_dict('records'):
        key = (row['Gene'], int(row['Position']))
        score = float(row['CAAS_score']) if row['CAAS_score'] not in ('NA', '') else None
        cur = out.get(key)
        if cur is None:
            out[key] = {'score': score, 'sides': [row['side']],
                        'p.emp': row['p.emp'], 'p.adj_bh': row['p.adj_bh'],
                        'caas_nominal': row['caas_nominal'], 'caas_adjusted': row['caas_adjusted'],
                        'evidence_status': row['evidence_status']}
        else:
            if score is not None and (cur['score'] is None or score > cur['score']):
                cur['score'] = score
            if row['side'] not in cur['sides']:
                cur['sides'].append(row['side'])
    return out


def build_position_characterization(positions, evidence, members):
    """One row per annotated or positively scored (Gene, Position)."""
    scores = position_scores_by_key(evidence)
    status_by_key = {(r.Gene, int(r.Position)): r.evaluation_status for r in positions.itertuples(index=False)}

    layers = {}
    for row in members.itertuples(index=False):
        key = (row.Gene, int(row.Position))
        layer = layers.setdefault(key, {'pfam_domain': set(), 'pfam_clan': set(),
                                        'ucr_region': None, 'variability': '',
                                        'fubar_call': '', 'fade_top': 0, 'fade_bottom': 0})
        detail = parse_detail(row.detail)
        if row.collection_id == 'pfam_domains':
            layer['pfam_domain'].add(row.term_label)
        elif row.collection_id == 'pfam_clans':
            layer['pfam_clan'].add(row.term_label)
        elif row.collection_id in ('ucr_core', 'ucr_flank'):
            region = detail.get('region', '')
            # 'core' wins over flank whenever both occur.
            if region == 'core' or layer['ucr_region'] in (None, '', 'flank_up', 'flank_down'):
                if layer['ucr_region'] != 'core':
                    layer['ucr_region'] = region
                    layer['variability'] = detail.get('variability', '')
        elif row.collection_id == 'fubar_tested':
            layer['fubar_call'] = detail.get('call', '')
        elif row.collection_id == 'fade_top_sig':
            layer['fade_top'] = 1
        elif row.collection_id == 'fade_bottom_sig':
            layer['fade_bottom'] = 1

    keys = sorted(set(layers) | set(scores))
    rows = []
    for gene, pos in keys:
        layer = layers.get((gene, pos), {'pfam_domain': set(), 'pfam_clan': set(),
                                        'ucr_region': None, 'variability': '',
                                        'fubar_call': '', 'fade_top': 0, 'fade_bottom': 0})
        sc = scores.get((gene, pos))
        rows.append({
            'Gene': gene, 'Position': pos,
            'evaluation_status': status_by_key.get((gene, pos), 'NA'),
            'pfam_domain': ';'.join(sorted(layer['pfam_domain'])),
            'pfam_clan': ';'.join(sorted(layer['pfam_clan'])),
            'ucr_region': layer['ucr_region'] or '',
            'position_variability': layer['variability'],
            'fubar_selection': layer['fubar_call'],
            'fade_top_sig': layer['fade_top'],
            'fade_bottom_sig': layer['fade_bottom'],
            'CAAS_score': 'NA' if sc is None or sc['score'] is None else repr(sc['score']),
            'sides': '' if sc is None else ';'.join(sorted(sc['sides'])),
            'p.emp': 'NA' if sc is None else sc['p.emp'],
            'p.adj_bh': 'NA' if sc is None else sc['p.adj_bh'],
            'caas_nominal': 'NA' if sc is None else sc['caas_nominal'],
            'caas_adjusted': 'NA' if sc is None else sc['caas_adjusted'],
        })
    return pd.DataFrame(rows, columns=['Gene', 'Position', 'evaluation_status', 'pfam_domain',
                                       'pfam_clan', 'ucr_region', 'position_variability',
                                       'fubar_selection', 'fade_top_sig', 'fade_bottom_sig',
                                       'CAAS_score', 'sides', 'p.emp', 'p.adj_bh',
                                       'caas_nominal', 'caas_adjusted'])


def build_collection_summary(collections, members, evidence, coverage):
    """Per collection/term: sizes and exact crosses with the ranking (no inference)."""
    scores = position_scores_by_key(evidence)
    n_annotable = 0
    for row in coverage.itertuples(index=False):
        if row.scope == 'universe' and row.metric == 'n_positions_annotable':
            n_annotable = int(float(row.value))
    rows = []
    for crow in collections.itertuples(index=False):
        if crow.status != 'included':
            continue
        sub = members[(members['collection_id'] == crow.collection_id) &
                      (members['term_id'] == crow.term_id)]
        keys = [(r.Gene, int(r.Position)) for r in sub.itertuples(index=False)]
        n_evaluated = sum(1 for k in keys if k in scores and scores[k]['score'] is not None)
        n_no_score = sum(1 for k in keys if k in scores and scores[k]['score'] is None)
        n_positive = sum(1 for k in keys if k in scores and scores[k]['score'] is not None)
        n_nominal = sum(1 for k in keys if k in scores and scores[k]['caas_nominal'] == 'TRUE')
        n_adjusted = sum(1 for k in keys if k in scores and scores[k]['caas_adjusted'] == 'TRUE')
        prop = (n_positive / n_evaluated) if n_evaluated else 'NA'
        rows.append({'collection_id': crow.collection_id, 'term_id': crow.term_id,
                     'term_type': crow.term_type, 'n_members': len(keys),
                     'n_genes': len({g for g, _ in keys}),
                     'n_members_evaluated': n_evaluated, 'n_members_without_canonical_score': n_no_score,
                     'n_members_positive': n_positive,
                     'n_members_nominal': n_nominal,
                     'n_members_adjusted': n_adjusted,
                     'prop_positive_of_evaluated': prop,
                     'positions_annotable_universe': n_annotable,
                     'denominator_definition': 'n_members_evaluated: miembros con score canónico positivo; no es el universo completo evaluado'})
    return pd.DataFrame(rows, columns=['collection_id', 'term_id', 'term_type', 'n_members',
                                       'n_genes', 'n_members_evaluated', 'n_members_without_canonical_score',
                                       'n_members_positive',
                                       'n_members_nominal', 'n_members_adjusted',
                                       'prop_positive_of_evaluated',
                                       'positions_annotable_universe', 'denominator_definition'])


def build_ranking_summary(evidence, n_cycles):
    """Score distribution and nominal/adjusted counts per ranking; no selection percentiles."""
    rows = []
    by_pos = {}
    for row in evidence.to_dict('records'):
        key = (row['Gene'], int(row['Position']))
        score = float(row['CAAS_score']) if row['CAAS_score'] not in ('NA', '') else None
        by_pos.setdefault(key, {'top': None, 'bottom': None, 'p.emp': row['p.emp'],
                                'caas_nominal': row['caas_nominal'], 'caas_adjusted': row['caas_adjusted']})
        if row['side'] in ('top', 'bottom') and score is not None:
            cur = by_pos[key][row['side']]
            if cur is None or score > cur:
                by_pos[key][row['side']] = score

    def summarize(ranking, values, keys):
        vals = [v for v in values if v is not None]
        nominal = sum(1 for k in keys if by_pos[k]['caas_nominal'] == 'TRUE')
        adjusted = sum(1 for k in keys if by_pos[k]['caas_adjusted'] == 'TRUE')
        p_unavailable = sum(1 for k in keys if by_pos[k]['p.emp'] in ('NA', ''))
        if vals:
            s = sorted(vals)
            q = lambda p: s[min(len(s) - 1, int(p * len(s)))]
            dist = {'min': repr(s[0]), 'q25': repr(q(0.25)), 'median': repr(q(0.5)),
                    'q75': repr(q(0.75)), 'max': repr(s[-1])}
        else:
            dist = {'min': 'NA', 'q25': 'NA', 'median': 'NA', 'q75': 'NA', 'max': 'NA'}
        rows.append({'ranking': ranking, 'n_positions': len(keys),
                     'n_scored': len(vals), **dist,
                     'n_nominal': nominal, 'n_adjusted': adjusted,
                     'n_p_unavailable': p_unavailable,
                     'N_cycles': 'NA' if n_cycles is None else str(int(n_cycles)),
                     'note': 'scores canónicos positivos disponibles; universo completo no disponible; nominal/adjusted son vistas'})

    # global: one unit per position (best side); top/bottom: directional rows.
    summarize('global', [max([v for v in (d['top'], d['bottom']) if v is not None], default=None)
                         for d in by_pos.values()], list(by_pos))
    for side in ('top', 'bottom'):
        keys = [k for k, d in by_pos.items() if d[side] is not None]
        summarize(side, [by_pos[k][side] for k in keys], keys)
    return pd.DataFrame(rows)


def build_universe_summary(positions, coverage):
    """Evaluation states of the position universe, with explicit denominators."""
    total = len(positions)
    rows = []
    if total:
        counts = positions['evaluation_status'].value_counts()
        for status, n in sorted(counts.items()):
            rows.append({'evaluation_status': status, 'n_positions': int(n),
                         'proportion': repr(n / total),
                         'denominator': total,
                         'denominator_definition': 'posiciones del registro (MAP + discovery + scores)'})
    else:
        rows.append({'evaluation_status': 'NA', 'n_positions': 0, 'proportion': 'NA',
                     'denominator': 0, 'denominator_definition': 'registro vacio'})
    for row in coverage.itertuples(index=False):
        if row.scope == 'universe':
            rows.append({'evaluation_status': f"coverage:{row.metric}", 'n_positions': row.value,
                         'proportion': 'NA', 'denominator': row.denominator,
                         'denominator_definition': row.denominator_definition})
    return pd.DataFrame(rows, columns=['evaluation_status', 'n_positions', 'proportion',
                                       'denominator', 'denominator_definition'])


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    missing = [f for f in REGISTRY_FILES
               if not os.path.exists(os.path.join(args.registry_dir, f))]
    if missing:
        raise SystemExit(f'registry incomplete, missing: {missing}')

    positions = read_tsv(os.path.join(args.registry_dir, 'positions.tsv'))
    evidence = read_tsv(os.path.join(args.registry_dir, 'position_evidence.tsv'))
    members = read_tsv(os.path.join(args.registry_dir, 'members.tsv'))
    collections = read_tsv(os.path.join(args.registry_dir, 'collections.tsv'))
    coverage = read_tsv(os.path.join(args.registry_dir, 'coverage.tsv'))
    with open(os.path.join(args.registry_dir, 'manifest.json')) as f:
        manifest = json.load(f)

    n_cycles = manifest.get('evidence', {}).get('null_cycles')
    n_cycles = None if n_cycles in ('NA', None) else float(n_cycles)

    pos_char = build_position_characterization(positions, evidence, members)
    coll_summary = build_collection_summary(collections, members, evidence, coverage)
    ranking = build_ranking_summary(evidence, n_cycles)
    universe = build_universe_summary(positions, coverage)

    write_tsv(pos_char, os.path.join(args.output_dir, 'position_characterization.tsv'))
    write_tsv(coll_summary, os.path.join(args.output_dir, 'collection_summary.tsv'))
    write_tsv(ranking, os.path.join(args.output_dir, 'ranking_summary.tsv'))
    write_tsv(universe, os.path.join(args.output_dir, 'universe_summary.tsv'))

    out_manifest = {
        'component': 'characterization',
        'contract_version': manifest.get('contract_version'),
        'inputs': {f: sha256_file(os.path.join(args.registry_dir, f)) for f in REGISTRY_FILES},
        'outputs': ['position_characterization.tsv', 'collection_summary.tsv',
                    'ranking_summary.tsv', 'universe_summary.tsv'],
        'thresholds': {'nominal': args.nominal_threshold, 'fdr': args.fdr_threshold,
                       'note': 'usados solo si la evidencia carece de flags canonicos'},
        'states': {
            'characterization': 'executed',
            'localization': 'not_executed',
            'global_association': 'not_executed',
            'proximity': 'not_executed',
        },
        'limits': ['sin inferencia de enriquecimiento',
                   'distribuciones de scores condicionadas a detecciones positivas; no son el universo completo evaluable',
                   'anotacion no implica significacion',
                   'p.emp compartido entre lados: no son dos pruebas direccionales',
                   'UCR/FUBAR/FADE no son validacion externa independiente'],
    }
    with open(os.path.join(args.output_dir, 'characterization_manifest.json'), 'w') as f:
        json.dump(out_manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    print(f"characterization: {len(pos_char)} characterized positions, "
          f"{len(coll_summary)} terms, {len(ranking)} rankings -> {args.output_dir}")


if __name__ == '__main__':
    main()
