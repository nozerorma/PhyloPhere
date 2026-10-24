#!/usr/bin/env python3
# posenrich_proximity.py — Proximity of prioritized positions to sites/regions (POSENRICH component 3).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Do the prioritized positions concentrate near sites or regions of interest?

WP10 P5 (docs/fases/WP10_posenrich_modular.md). Protein and genomic proximity are
separate contracts: protein distance is measured in reference residues
(hg38_aa_pos), never in MSA columns with gaps; genomic distance is measured in
nucleotides and never mixed into the protein result.

Protein proximity: distance in reference residues from each evaluated position to
the nearest member of the same protein of the term; 0 inside a region or on the
site itself. Proteins without a comparable reference are reported as coverage
lack, never as infinite distance. Query positions carry canonical positive
scores; no p-value filtering is applied. This is not the full evaluable universe.

Genomic proximity: distance in nucleotides on the same chromosome to the nearest
member of the term; 0 on coincidence. Window membership (1 Mb bins) stays in the
registry members and is not re-measured here.

The distance control is not assumed to be the AUC control: the first
implementation is descriptive with an explicit state; radii are fixed before any
contrast and reported as descriptive.

Outputs (in --output_dir):
  proximity_protein_distances.tsv   one row per evaluated position x term (genes with sites)
  proximity_protein_summary.tsv     per collection/term: distribution summary + fixed radii
  proximity_genomic_distances.tsv   one row per evaluated position x site term (same chromosome)
  proximity_genomic_summary.tsv     per collection/term summary
  proximity_manifest.json           parameters, states and limits
"""

import argparse
import csv
import json
import os

import pandas as pd
from posenrich_localization import position_ranking_scores

REGISTRY_FILES = ['positions.tsv', 'position_evidence.tsv', 'members.tsv',
                  'collections.tsv', 'coverage.tsv', 'manifest.json']

# Collections measured as protein sites/regions; coverage and genomic-window
# collections are not protein proximity targets.
PROTEIN_COLLECTIONS = ('pfam_domains', 'pfam_clans', 'ucr_core', 'ucr_flank',
                       'fubar_positive', 'fubar_purifying', 'fade_top_sig',
                       'fade_bottom_sig', 'custom')
REGION_TYPES = {'pfam_domains', 'pfam_clans', 'ucr_core', 'ucr_flank'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--registry_dir', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--collections', default=None,
                   help='Optional comma-separated collection_id filter (protein proximity)')
    p.add_argument('--radios', default='0,5,10',
                   help='Fixed descriptive radii in reference residues, fixed before any contrast')
    return p.parse_args()


def read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')


def parse_coord(coord):
    """'chr1:1000005' -> ('chr1', 1000005); None when not a genomic coordinate."""
    if coord in ('NA', '', None):
        return None
    parts = str(coord).split(':')
    if len(parts) != 2 or not parts[0].startswith('chr'):
        return None
    try:
        return parts[0], int(parts[1])
    except ValueError:
        return None


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    missing = [f for f in REGISTRY_FILES
               if not os.path.exists(os.path.join(args.registry_dir, f))]
    if missing:
        raise SystemExit(f'registry incomplete, missing: {missing}')
    radios = sorted({int(r) for r in args.radios.split(',') if r.strip() != ''})
    if any(r < 0 for r in radios):
        raise SystemExit('--radios must be >= 0')

    positions = read_tsv(os.path.join(args.registry_dir, 'positions.tsv'))
    evidence = read_tsv(os.path.join(args.registry_dir, 'position_evidence.tsv'))
    members = read_tsv(os.path.join(args.registry_dir, 'members.tsv'))
    collections = read_tsv(os.path.join(args.registry_dir, 'collections.tsv'))
    with open(os.path.join(args.registry_dir, 'manifest.json')) as f:
        manifest = json.load(f)

    # Coordinates by (gene, pos); evaluated positions carry their score.
    coords = {}
    for row in positions.to_dict('records'):
        coords[(row['Gene'], int(row['Position']))] = {
            'residue': None if row['hg38_aa_pos'] in ('NA', '') else int(row['hg38_aa_pos']),
            'nt': parse_coord(row['hg38_nt_coord']),
            'in_map': row['in_map'] == '1',
        }
    evaluated = sorted({(r['Gene'], int(r['Position'])) for r in evidence.to_dict('records')
                        if r['CAAS_score'] not in ('NA', '')})
    scores = position_ranking_scores(evidence, 'caas_score')

    site_terms = collections[(collections['status'] == 'included') &
                             (collections['collection_id'].isin(PROTEIN_COLLECTIONS))]
    if args.collections:
        wanted = set(args.collections.split(','))
        site_terms = site_terms[site_terms['collection_id'].isin(wanted)]
    region_type = dict(zip(collections['collection_id'] + '|' + collections['term_id'],
                           collections['term_type']))

    # term -> gene -> member positions
    term_members = {}
    for row in members.to_dict('records'):
        cid = row['collection_id']
        if cid not in PROTEIN_COLLECTIONS:
            continue
        if cid not in set(site_terms['collection_id']):
            continue
        term_members.setdefault((cid, row['term_id']), {}).setdefault(
            row['Gene'], set()).add(int(row['Position']))

    protein_rows = []
    genomic_rows = []
    for (cid, term), by_gene in sorted(term_members.items()):
        is_region = region_type.get(f'{cid}|{term}') == 'region'
        for gene in sorted(by_gene):
            site_positions = by_gene[gene]
            site_residues = [coords[(gene, p)]['residue'] for p in sorted(site_positions)
                             if (gene, p) in coords]
            site_residues = [r for r in site_residues if r is not None]
            site_nt = [coords[(gene, p)]['nt'] for p in sorted(site_positions)
                       if (gene, p) in coords]
            site_nt = [c for c in site_nt if c is not None]
            for (g, pos) in evaluated:
                if g != gene:
                    continue
                info = coords.get((g, pos))
                if info is None:
                    continue
                score = scores.get((g, pos), 'NA')
                # protein distance
                if info['residue'] is None:
                    protein_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                         'Position': pos, 'hg38_aa_pos': 'NA',
                                         'distance_residues': 'NA', 'in_region': 'NA',
                                         'CAAS_score': score,
                                         'status': 'no_reference_residue'})
                elif not site_residues:
                    protein_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                         'Position': pos, 'hg38_aa_pos': info['residue'],
                                         'distance_residues': 'NA', 'in_region': 'NA',
                                         'CAAS_score': score,
                                         'status': 'no_sites_with_reference'})
                else:
                    in_region = 1 if (is_region and pos in site_positions) else 0
                    if pos in site_positions:
                        distance = 0
                    else:
                        distance = min(abs(info['residue'] - r) for r in site_residues)
                    protein_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                         'Position': pos, 'hg38_aa_pos': info['residue'],
                                         'distance_residues': distance,
                                         'in_region': in_region,
                                         'CAAS_score': score, 'status': 'measured'})
                # genomic distance (separate contract)
                if info['nt'] is None:
                    genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                         'Position': pos, 'hg38_nt_coord': 'NA',
                                         'distance_nt': 'NA', 'CAAS_score': score,
                                         'status': 'no_genomic_coordinate'})
                elif not site_nt:
                    genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                         'Position': pos, 'hg38_nt_coord': f'{info["nt"][0]}:{info["nt"][1]}',
                                         'distance_nt': 'NA', 'CAAS_score': score,
                                         'status': 'no_sites_with_genomic_coordinate'})
                else:
                    same_chrom = [c for c in site_nt if c[0] == info['nt'][0]]
                    if not same_chrom:
                        genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                             'Position': pos,
                                             'hg38_nt_coord': f'{info["nt"][0]}:{info["nt"][1]}',
                                             'distance_nt': 'NA', 'CAAS_score': score,
                                             'status': 'no_sites_same_chromosome'})
                    else:
                        distance = 0 if pos in site_positions else min(
                            abs(info['nt'][1] - c[1]) for c in same_chrom)
                        genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                             'Position': pos,
                                             'hg38_nt_coord': f'{info["nt"][0]}:{info["nt"][1]}',
                                             'distance_nt': distance,
                                             'CAAS_score': score, 'status': 'measured'})

    protein_df = pd.DataFrame(protein_rows, columns=[
        'collection_id', 'term_id', 'Gene', 'Position', 'hg38_aa_pos', 'distance_residues',
        'in_region', 'CAAS_score', 'status'])
    genomic_df = pd.DataFrame(genomic_rows, columns=[
        'collection_id', 'term_id', 'Gene', 'Position', 'hg38_nt_coord', 'distance_nt',
        'CAAS_score', 'status'])

    def summarize(df, dist_col, summary_path):
        rows = []
        for (cid, term), sub in df.groupby(['collection_id', 'term_id']):
            measured = sub[sub['status'] == 'measured']
            n_measured = len(measured)
            n_evaluated = int((sub['status'] != 'measured').sum()) + n_measured
            distances = pd.to_numeric(measured[dist_col], errors='coerce')
            row = {'collection_id': cid, 'term_id': term,
                   'n_positions_evaluated': n_evaluated,
                   'n_positions_measured': n_measured,
                   'coverage_proportion': repr(n_measured / n_evaluated) if n_evaluated else 'NA',
                   'median_distance': repr(float(distances.median())) if n_measured else 'NA',
                   'mean_distance': repr(float(distances.mean())) if n_measured else 'NA',
                   'n_genes_with_sites': int(sub['Gene'].nunique()),
                   'control_status': 'not_validated',
                   'status': 'descriptive'}
            for radius in radios:
                prop = float((distances <= radius).mean()) if n_measured else 'NA'
                row[f'prop_within_{radius}'] = repr(prop)
            rows.append(row)
        cols = ['collection_id', 'term_id', 'n_positions_evaluated', 'n_positions_measured',
                'coverage_proportion', 'median_distance', 'mean_distance', 'n_genes_with_sites'] \
            + [f'prop_within_{r}' for r in radios] + ['control_status', 'status']
        out = pd.DataFrame(rows, columns=cols)
        write_tsv(out, summary_path)
        return out

    summarize(protein_df, 'distance_residues',
              os.path.join(args.output_dir, 'proximity_protein_summary.tsv'))
    summarize(genomic_df, 'distance_nt',
              os.path.join(args.output_dir, 'proximity_genomic_summary.tsv'))
    write_tsv(protein_df, os.path.join(args.output_dir, 'proximity_protein_distances.tsv'))
    write_tsv(genomic_df, os.path.join(args.output_dir, 'proximity_genomic_distances.tsv'))

    out_manifest = {
        'component': 'proximity',
        'ranking_population': 'canonical_positive_detections_only',
        'contract_version': manifest.get('contract_version'),
        'parameters': {'collections': args.collections, 'radios': radios},
        'contracts': {'protein': 'distancia en residuos de referencia (hg38_aa_pos); 0 dentro de region',
                      'genomic': 'distancia en nt en el mismo cromosoma; 0 en coincidencia'},
        'states': {'proximity': 'executed_descriptive' if len(protein_df) else 'no_terms',
                   'characterization': 'not_executed (componente separado)',
                   'localization': 'not_executed (componente separado)',
                   'global_association': 'not_executed (componente separado)'},
        'control': {'status': 'not_validated',
                    'reason': 'el control de AUC no se asume valido para distancias; '
                              'no se infiere p en esta primera implementacion'},
        'limits': ['radios descriptivos fijados antes del contraste',
                   'proteinas sin referencia comparable se informan como falta de cobertura, no distancia infinita',
                   'proximidad genómica y proteica no se mezclan',
                   'la distancia no se mide en columnas MSA con gaps',
                   'consulta condicionada a detecciones positivas; no es el universo completo evaluado'],
    }
    with open(os.path.join(args.output_dir, 'proximity_manifest.json'), 'w') as f:
        json.dump(out_manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    print(f"proximity: {len(protein_df)} protein distance rows, {len(genomic_df)} genomic rows "
          f"-> {args.output_dir}")


if __name__ == '__main__':
    main()
