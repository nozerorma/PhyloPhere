#!/usr/bin/env python3
# posenrich_global.py — Global association of position sets with the ranking (POSENRICH component 2b).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Does a position set rank higher than the compatible global background?

WP10 P4 (docs/fases/WP10_posenrich_modular.md). This is the global association
question, kept separate from within-protein localization: a set may accumulate
signal simply because certain proteins accumulate signal. Nothing here reuses the
localization p-value.

Estimand: AUC of the set members against every other evaluable position of the
universe for the ranking (no hit selection before the AUC). Composition by gene
and leave-one-gene-out sensitivity are always reported.

Calibration: a global p needs the statistic reproduced per null cycle. The
positional null is kept as counters; marginal counters do not reconstruct scores
or co-occurrences per cycle, and the stored max-over-sides p.emp does not support
a top/bottom null. Without a reproducible per-cycle statistic the result stays
descriptive with an explicit state; no substitute p is fabricated.

Outputs (in --output_dir):
  global_results.tsv   long table: one row per (ranking, collection, term)
  global_by_gene.tsv   composition per gene and leave-one-gene-out sensitivity
  global_manifest.json parameters, states and limits
"""

import argparse
import csv
import json
import os

import numpy as np
import pandas as pd

REGISTRY_FILES = ['positions.tsv', 'position_evidence.tsv', 'members.tsv',
                  'collections.tsv', 'coverage.tsv', 'manifest.json']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--registry_dir', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--ranking', default='caas_score',
                   help="caas_score (best-side continuous score, principal) or "
                        "minus_log10_pemp (explicit sensitivity over the shared p.emp)")
    p.add_argument('--collections', default=None,
                   help='Optional comma-separated collection_id filter')
    return p.parse_args()


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


def position_ranking_scores(evidence, ranking):
    """(gene, pos) -> score; caas_score keeps the best side per position."""
    import math
    out = {}
    for row in evidence.to_dict('records'):
        key = (row['Gene'], int(row['Position']))
        if ranking == 'minus_log10_pemp':
            p = row['p.emp']
            if p in ('NA', '') or float(p) <= 0:
                continue
            score = -math.log10(float(p))
        else:
            v = row['CAAS_score']
            if v in ('NA', ''):
                continue
            score = float(v)
        cur = out.get(key)
        if cur is None or score > cur:
            out[key] = score
    return out


def average_ranks(scores):
    order = np.argsort(scores, kind='mergesort')
    ranks = np.empty(len(scores), dtype=float)
    sorted_scores = scores[order]
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def auc(scores, member_mask):
    """AUC with half credit to ties; None without both groups."""
    n1 = int(member_mask.sum())
    n0 = len(member_mask) - n1
    if n1 == 0 or n0 == 0:
        return None
    ranks = average_ranks(scores)
    return (float(ranks[member_mask == 1].sum()) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    missing = [f for f in REGISTRY_FILES
               if not os.path.exists(os.path.join(args.registry_dir, f))]
    if missing:
        raise SystemExit(f'registry incomplete, missing: {missing}')
    if args.ranking not in ('caas_score', 'minus_log10_pemp'):
        raise SystemExit(f'unknown ranking: {args.ranking}')

    evidence = read_tsv(os.path.join(args.registry_dir, 'position_evidence.tsv'))
    members = read_tsv(os.path.join(args.registry_dir, 'members.tsv'))
    collections = read_tsv(os.path.join(args.registry_dir, 'collections.tsv'))
    with open(os.path.join(args.registry_dir, 'manifest.json')) as f:
        manifest = json.load(f)

    scores = position_ranking_scores(evidence, args.ranking)
    keys = sorted(scores)
    score_arr = np.array([scores[k] for k in keys], dtype=float)
    member_keys = set(zip(members['Gene'], members['Position'].astype(int)))
    mask_all = np.array([1 if k in member_keys else 0 for k in keys], dtype=float)

    eligible = collections[collections['status'] == 'included']
    if args.collections:
        wanted = set(args.collections.split(','))
        eligible = eligible[eligible['collection_id'].isin(wanted)]

    result_rows = []
    gene_rows = []
    for row in eligible.to_dict('records'):
        cid, term = row['collection_id'], row['term_id']
        sub = members[(members['collection_id'] == cid) & (members['term_id'] == term)]
        term_keys = [(r['Gene'], int(r['Position'])) for r in sub.to_dict('records')]
        # Members outside the evaluable background are audited, not dropped silently.
        n_members_total = len(term_keys)
        member_set = set(term_keys)
        mask = np.array([1 if k in member_set else 0 for k in keys], dtype=float)
        n_members = int(mask.sum())
        auc_obs = auc(score_arr, mask)
        genes = sorted({g for g, _ in term_keys})
        n_genes = len(genes)
        if auc_obs is None:
            result_rows.append({
                'analysis': 'global_association', 'scope': 'global', 'ranking': args.ranking,
                'collection_id': cid, 'term_id': term, 'effect_statistic': 'AUC',
                'auc': 'NA', 'n_members': n_members, 'n_members_total': n_members_total,
                'n_background': len(keys) - n_members, 'n_genes': n_genes,
                'loo_auc_min': 'NA', 'loo_auc_max': 'NA',
                'calibration_status': 'not_available', 'p': 'NA', 'status': 'not_identifiable',
                'reason': 'el conjunto no tiene miembros evaluables o no tiene fondo',
            })
            continue
        # Leave-one-gene-out sensitivity: composition and influence of each protein.
        loo = []
        for gene in genes:
            keep = [i for i, k in enumerate(keys) if k[0] != gene]
            if not keep:
                continue
            loo_auc = auc(score_arr[keep], mask[keep])
            if loo_auc is not None:
                loo.append((gene, loo_auc))
                gene_rows.append({'collection_id': cid, 'term_id': term, 'ranking': args.ranking,
                                  'Gene': gene, 'n_members': int(sum(1 for g, _ in term_keys if g == gene)),
                                  'auc_without_gene': repr(loo_auc),
                                  'delta_auc': repr(auc_obs - loo_auc)})
        loo_values = [v for _, v in loo]
        result_rows.append({
            'analysis': 'global_association', 'scope': 'global', 'ranking': args.ranking,
            'collection_id': cid, 'term_id': term, 'effect_statistic': 'AUC',
            'auc': repr(auc_obs), 'n_members': n_members, 'n_members_total': n_members_total,
            'n_background': len(keys) - n_members, 'n_genes': n_genes,
            'loo_auc_min': repr(min(loo_values)) if loo_values else 'NA',
            'loo_auc_max': repr(max(loo_values)) if loo_values else 'NA',
            'calibration_status': 'not_available',
            'p': 'NA',
            'status': 'descriptive',
            'reason': 'sin calibracion global reproducible por ciclo; los contadores marginales '
                      'no reconstruyen scores ni coocurrencias; no se sustituye por el p de localizacion',
        })

    write_tsv(pd.DataFrame(result_rows), os.path.join(args.output_dir, 'global_results.tsv'))
    write_tsv(pd.DataFrame(gene_rows), os.path.join(args.output_dir, 'global_by_gene.tsv'))

    out_manifest = {
        'component': 'global_association',
        'scope': 'global',
        'contract_version': manifest.get('contract_version'),
        'parameters': {'ranking': args.ranking, 'collections': args.collections},
        'calibration': {'status': 'not_available',
                        'reason': 'el estadistico por ciclo no es reproducible desde contadores marginales',
                        'p_localization_reused': False},
        'states': {'global_association': 'executed_descriptive' if result_rows else 'no_terms',
                   'localization': 'not_executed (componente separado)',
                   'characterization': 'not_executed (componente separado)',
                   'proximity': 'not_executed (componente separado)'},
        'limits': ['asociacion global, no contraste dentro de proteina',
                   'un conjunto puede acumular señal solo porque ciertas proteinas acumulan señal',
                   'sin p: la calibracion por fenotipos requiere el estadistico por ciclo'],
    }
    with open(os.path.join(args.output_dir, 'global_manifest.json'), 'w') as f:
        json.dump(out_manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    print(f"global: {len(result_rows)} terms, ranking={args.ranking} -> {args.output_dir}")


if __name__ == '__main__':
    main()
