#!/usr/bin/env python3
# posenrich_localization.py — Within-protein localization contrast (POSENRICH component 2, priority).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Do the members of a position set rank higher than the other eligible positions
of the same proteins?

WP10 P3 (docs/fases/WP10_posenrich_modular.md). This is the within-protein
localization contrast: association global (which may just accumulate signal in
certain proteins) is a separate component (posenrich_global.py).

Estimand (fixed before the engine, see docs/POSENRICH_MODULAR_CONTRATO.md):
  * Unit: one biological position (Gene, CAAS 0-based Position); never multiplied
    by side, hypothesis or evidence rows. The pooled p.emp is shared between sides.
  * Ranking: --ranking caas_score (default; best-side continuous CAAS score over
    evaluable positions, no hit selection before the AUC) or minus_log10_pemp
    (explicit sensitivity over the shared p.emp). top/bottom select the
    directional score.
  * Within-gene effect: AUC of members vs eligible non-members of the same
    protein, half credit to ties (average ranks). Only genes with both groups
    inform.
  * Aggregation: --aggregation mean (default principal; equal weight per gene)
    or pair_weighted (sensitivity). Both columns are always reported.
  * Control: within-protein label redistribution. Membership labels are permuted
    inside each protein, preserving protein, score profile and coverage; no
    uniform shuffle across proteins, no circular shifts. p Monte Carlo unilateral
    for greater AUC: inclusive exceedances, (k+1)/(B_valid+1), fixed seed.

Outputs (in --output_dir):
  localization_results.tsv   long table: one row per (ranking, collection, term)
  localization_by_gene.tsv   per-gene AUC detail (informative genes only)
  localization_manifest.json parameters, seed, B, states and limits
"""

import argparse
import csv
import json
import math
import os

import numpy as np
import pandas as pd

REGISTRY_FILES = ['positions.tsv', 'position_evidence.tsv', 'members.tsv',
                  'collections.tsv', 'coverage.tsv', 'manifest.json']

RANKINGS = ('caas_score', 'minus_log10_pemp')
# caas_score is the continuous CAAS score: 'global' keeps the best side per position,
# 'top'/'bottom' keep the directional score of their side.
SCORE_SCOPES = ('global', 'top', 'bottom')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--registry_dir', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--ranking', default='caas_score', choices=RANKINGS,
                   help='caas_score: continuous best-side score (principal); '
                        'minus_log10_pemp: explicit sensitivity over the shared p.emp')
    p.add_argument('--aggregation', default='mean', choices=('mean', 'pair_weighted'),
                   help='mean: equal weight per informative gene (principal); '
                        'pair_weighted: sensitivity by n_members*n_background pairs')
    p.add_argument('--n_permutations', type=int, default=999)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--collections', default=None,
                   help='Optional comma-separated collection_id filter')
    return p.parse_args()


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


def average_ranks(scores):
    """Average ranks (1-based) with half credit to ties."""
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


def auc_from_ranks(ranks, member_mask):
    """AUC = (sum of member ranks - n1(n1+1)/2) / (n1*n0); None without both groups."""
    n1 = int(member_mask.sum())
    n0 = len(member_mask) - n1
    if n1 == 0 or n0 == 0:
        return None
    return (float(ranks[member_mask == 1].sum()) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def position_ranking_scores(evidence, ranking, scope='global'):
    """(gene, pos) -> score for the ranking, or None when not evaluable.

    caas_score/global: best-side score (one unit per position). caas_score/top|bottom:
    the directional score of that side. minus_log10_pemp: -log10 of the shared
    max-over-sides p.emp; positions without p.emp are not evaluable for it.
    """
    out = {}
    for row in evidence.to_dict('records'):
        key = (row['Gene'], int(row['Position']))
        if ranking == 'minus_log10_pemp':
            p = row['p.emp']
            if p in ('NA', '') or float(p) <= 0:
                continue
            score = -math.log10(float(p))
        else:
            if scope in ('top', 'bottom') and row['side'] != scope:
                continue
            v = row['CAAS_score']
            if v in ('NA', ''):
                continue
            score = float(v)
        cur = out.get(key)
        if cur is None or score > cur:
            out[key] = score
    return out


def build_gene_table(members, collections, scores, collection_filter):
    """(collection_id, term_id) -> list of per-gene (scores array, member mask, n_members_total)."""
    eligible = collections[(collections['informative_within_protein'] == '1') &
                           (collections['status'] == 'included')]
    if collection_filter:
        wanted = set(collection_filter.split(','))
        eligible = eligible[eligible['collection_id'].isin(wanted)]
    members = members[members['collection_id'].isin(set(eligible['collection_id']))]
    member_keys = set(zip(members['Gene'], members['Position'].astype(int)))
    # evaluable positions per gene: those with a score in this ranking
    by_gene = {}
    for (gene, pos), score in scores.items():
        by_gene.setdefault(gene, {})[pos] = score
    terms = {}
    for row in eligible.to_dict('records'):
        cid, term = row['collection_id'], row['term_id']
        genes = {}
        sub = members[(members['collection_id'] == cid) & (members['term_id'] == term)]
        member_by_gene = {}
        for r in sub.to_dict('records'):
            member_by_gene.setdefault(r['Gene'], set()).add(int(r['Position']))
        n_members_total = sum(len(v) for v in member_by_gene.values())
        for gene, positions in by_gene.items():
            member_positions = member_by_gene.get(gene, set())
            # Members outside the evaluable set are audited, not dropped silently.
            n_members_evaluable = sum(1 for p in member_positions if p in positions)
            genes[gene] = {'positions': positions, 'member_positions': member_positions,
                           'n_members_evaluable': n_members_evaluable}
        terms[(cid, term)] = {'genes': genes, 'n_members_total': n_members_total,
                              'description': row['description']}
    return terms


def run_term(term_data, aggregation, n_permutations, rng):
    """AUC per informative gene, aggregation, and the within-protein permutation control."""
    per_gene = []
    prepared = []
    for gene in sorted(term_data['genes']):
        info = term_data['genes'][gene]
        positions = info['positions']
        member_positions = info['member_positions']
        pos_list = sorted(positions)
        score_arr = np.array([positions[p] for p in pos_list], dtype=float)
        mask = np.array([1 if p in member_positions else 0 for p in pos_list], dtype=float)
        auc = auc_from_ranks(average_ranks(score_arr), mask)
        if auc is None:
            continue  # gene without both groups: not informative
        ranks = average_ranks(score_arr)
        n1 = int(mask.sum())
        prepared.append((gene, score_arr, mask, ranks, n1, auc))
        per_gene.append({'Gene': gene, 'auc': auc, 'n_members': n1,
                         'n_background': len(pos_list) - n1,
                         'n_positions_evaluable': len(pos_list)})
    if not prepared:
        return None, per_gene, 0

    aucs = np.array([p[5] for p in prepared])
    weights = np.array([p[4] * (len(p[1]) - p[4]) for p in prepared], dtype=float)

    def aggregate(values):
        if aggregation == 'pair_weighted':
            return float((values * weights).sum() / weights.sum())
        return float(values.mean())

    auc_obs = aggregate(aucs)
    auc_mean = float(aucs.mean())
    auc_pair = float((aucs * weights).sum() / weights.sum())

    # Within-protein label redistribution: per gene, uniform permutation of the
    # membership labels (fixed n1 per gene), preserving scores and coverage.
    perm_aucs = np.empty(n_permutations, dtype=float)
    b_valid = 0
    for b in range(n_permutations):
        values = np.empty(len(prepared), dtype=float)
        for gi, (gene, score_arr, mask, ranks, n1, _auc) in enumerate(prepared):
            n = len(score_arr)
            order = rng.permutation(n)
            perm_mask = np.zeros(n)
            perm_mask[order[:n1]] = 1.0
            values[gi] = auc_from_ranks(ranks, perm_mask)
        perm_aucs[b] = aggregate(values)
        b_valid += 1
    k = int((perm_aucs >= auc_obs - 1e-12).sum())
    p_mc = (k + 1) / (b_valid + 1)
    return {'auc': auc_obs, 'auc_mean': auc_mean, 'auc_pair_weighted': auc_pair,
            'auc_control_mean': float(perm_aucs.mean()),
            'p_mc': p_mc, 'k_exceedances': k, 'B_valid': b_valid}, per_gene, len(prepared)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    missing = [f for f in REGISTRY_FILES
               if not os.path.exists(os.path.join(args.registry_dir, f))]
    if missing:
        raise SystemExit(f'registry incomplete, missing: {missing}')
    if args.n_permutations < 1:
        raise SystemExit('--n_permutations must be >= 1')

    evidence = read_tsv(os.path.join(args.registry_dir, 'position_evidence.tsv'))
    members = read_tsv(os.path.join(args.registry_dir, 'members.tsv'))
    collections = read_tsv(os.path.join(args.registry_dir, 'collections.tsv'))
    with open(os.path.join(args.registry_dir, 'manifest.json')) as f:
        manifest = json.load(f)
    n_cycles = manifest.get('evidence', {}).get('null_cycles')

    scores = position_ranking_scores(evidence, args.ranking)
    terms = build_gene_table(members, collections, scores, args.collections)
    rng = np.random.default_rng(args.seed)

    result_rows = []
    gene_rows = []
    n_not_identifiable = 0
    for (cid, term) in sorted(terms):
        term_data = terms[(cid, term)]
        summary, per_gene, n_informative = run_term(term_data, args.aggregation,
                                                    args.n_permutations, rng)
        n_genes_evaluated = len(term_data['genes'])
        if summary is None:
            n_not_identifiable += 1
            result_rows.append({
                'analysis': 'localization', 'scope': 'within_protein', 'ranking': args.ranking,
                'collection_id': cid, 'term_id': term, 'effect_statistic': 'AUC',
                'auc': 'NA', 'auc_minus_half': 'NA', 'auc_control_mean': 'NA',
                'auc_excess_vs_control': 'NA', 'p_mc': 'NA',
                'n_genes_informative': 0, 'n_genes_evaluated': n_genes_evaluated,
                'n_members': term_data['n_members_total'], 'n_background_positions': 'NA',
                'n_permutations': args.n_permutations, 'B_valid': 0, 'seed': args.seed,
                'aggregation': args.aggregation, 'coverage_status': 'no_within_protein_contrast',
                'status': 'not_identifiable',
                'reason': 'ningun gen tiene miembros y fondo evaluables con contraste',
            })
            continue
        n_background = sum(g['n_background'] for g in per_gene)
        result_rows.append({
            'analysis': 'localization', 'scope': 'within_protein', 'ranking': args.ranking,
            'collection_id': cid, 'term_id': term, 'effect_statistic': 'AUC',
            'auc': repr(summary['auc']), 'auc_mean': repr(summary['auc_mean']),
            'auc_pair_weighted': repr(summary['auc_pair_weighted']),
            'auc_minus_half': repr(summary['auc'] - 0.5),
            'auc_control_mean': repr(summary['auc_control_mean']),
            'auc_excess_vs_control': repr(summary['auc'] - summary['auc_control_mean']),
            'p_mc': repr(summary['p_mc']),
            'n_genes_informative': n_informative, 'n_genes_evaluated': n_genes_evaluated,
            'n_members': term_data['n_members_total'], 'n_background_positions': n_background,
            'n_permutations': args.n_permutations, 'B_valid': summary['B_valid'],
            'seed': args.seed, 'aggregation': args.aggregation,
            'coverage_status': 'evaluated',
            'status': 'executed',
            'reason': '',
        })
        for g in per_gene:
            gene_rows.append({'collection_id': cid, 'term_id': term, 'ranking': args.ranking,
                              'Gene': g['Gene'], 'auc': repr(g['auc']),
                              'n_members': g['n_members'], 'n_background': g['n_background'],
                              'n_positions_evaluable': g['n_positions_evaluable']})

    results = pd.DataFrame(result_rows)
    write_tsv(results, os.path.join(args.output_dir, 'localization_results.tsv'))
    write_tsv(pd.DataFrame(gene_rows), os.path.join(args.output_dir, 'localization_by_gene.tsv'))

    out_manifest = {
        'component': 'localization',
        'scope': 'within_protein',
        'contract_version': manifest.get('contract_version'),
        'parameters': {'ranking': args.ranking, 'aggregation': args.aggregation,
                       'n_permutations': args.n_permutations, 'seed': args.seed,
                       'collections': args.collections},
        'control': {'type': 'within_protein_label_redistribution',
                    'p_mc': 'unilateral greater AUC, inclusive exceedances, (k+1)/(B_valid+1)',
                    'preserves': ['protein', 'score profile', 'coverage', 'members per gene']},
        'states': {'localization': 'executed' if len(result_rows) else 'no_terms',
                    'characterization': 'not_executed (separate component)',
                    'global_association': 'not_executed (separate component)',
                    'proximity': 'not_executed (separate component)'},
        'terms_not_identifiable': n_not_identifiable,
        'limits': ['el control condicionado no valida la asociacion global',
                   'AUC 0.5 no es necesariamente el centro del control condicionado',
                   'no selecciona hits antes del AUC',
                   'p.emp compartido entre lados; el score direccional es aparte'],
    }
    with open(os.path.join(args.output_dir, 'localization_manifest.json'), 'w') as f:
        json.dump(out_manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    print(f"localization: {len(result_rows)} terms ({n_not_identifiable} not identifiable), "
          f"ranking={args.ranking}, B={args.n_permutations}, seed={args.seed} -> {args.output_dir}")


if __name__ == '__main__':
    main()
