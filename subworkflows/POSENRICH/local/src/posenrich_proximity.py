#!/usr/bin/env python3
# posenrich_proximity.py: proximity of prioritized positions to sites and regions (POSENRICH).
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""Do the prioritized positions concentrate near sites or regions of interest?

Protein and genomic proximity are
separate contracts: protein distance is measured in reference residues
(hg38_aa_pos), never in MSA columns with gaps; genomic distance is measured in
nucleotides and never mixed into the protein result.

Protein proximity: distance in reference residues from each evaluated position to
the nearest member of the same protein of the term; 0 inside a region or on the
site itself. Proteins without a comparable reference are reported as coverage
lack, never as infinite distance. Query positions carry canonical positive
scores; no p-value filtering is applied. This is not the full evaluable universe.

Genomic proximity: distance in nucleotides to the nearest member of the term
within the same gene and chromosome; 0 on coincidence. Window membership (1 Mb bins) stays in the
registry members and is not re-measured here.

No p-value is inferred: results are descriptive and radii are fixed before any contrast.

Protein hits list the detections that coincide with a member of a term (relation overlap, by identity of the
position) or lie within the largest radius of one (relation near), with the member they are closest to; the
by-gene table counts them per term and gene.

Outputs (in --output_dir):
  proximity_protein_distances.tsv   one row per evaluated position x term (genes with sites)
  proximity_protein_summary.tsv     per collection/term: distribution summary + fixed radii
  proximity_protein_hits.tsv        overlaps and near detections, one row per detection x term
  proximity_protein_by_gene.tsv     per term and gene: members, detections, overlaps, near detections
  proximity_genomic_distances.tsv   one row per evaluated position x site term (same gene/chromosome)
  proximity_genomic_summary.tsv     per collection/term summary
  proximity_manifest.json           parameters, states and limits
"""

import argparse
import csv
import json
import os
from array import array
from bisect import bisect_left
from contextlib import ExitStack

import numpy as np

import posenrich_members as pm

REGISTRY_FILES = ['positions.tsv', 'position_evidence.tsv', 'regions.tsv', 'sites.tsv',
                  'collections.tsv', 'coverage.tsv', 'manifest.json']

# Collections measured as protein sites/regions; coverage and genomic-window
# collections are not protein proximity targets.
PROTEIN_COLLECTIONS = ('pfam_domains', 'pfam_clans', 'ucr_core', 'ucr_flank',
                       'fubar_positive', 'fubar_purifying', 'fade_top_sig',
                       'fade_bottom_sig', 'custom')
REGION_TYPES = {'pfam_domains', 'pfam_clans', 'ucr_core', 'ucr_flank'}
HIT_COLUMNS = ['collection_id', 'term_id', 'Gene', 'Position', 'hg38_aa_pos', 'relation', 'distance_residues',
               'nearest_site_position', 'nearest_site_residue', 'in_region', 'CAAS_score']
BY_GENE_COLUMNS = ['collection_id', 'term_id', 'Gene', 'n_term_positions_in_gene', 'n_detections_measured',
                   'n_overlaps', 'n_near', 'min_distance_residues']


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
    with open(path, newline='') as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def nearest_site(value, sorted_sites):
    """(distance, site) of the nearest site in logarithmic time, including gaps and disjoint regions."""
    index = bisect_left(sorted_sites, value)
    return min((abs(value - sorted_sites[i]), sorted_sites[i]) for i in (index - 1, index)
               if 0 <= i < len(sorted_sites))


def nearest_distance(value, sorted_sites):
    return nearest_site(value, sorted_sites)[0]


class DistanceOutput:
    """Stream detail rows; retain compact distances for one term's exact median."""

    def __init__(self, stack, output_dir, kind, coord_col, dist_col, radios):
        self.dist_col = dist_col
        self.radios = radios
        columns = ['collection_id', 'term_id', 'Gene', 'Position', coord_col, dist_col]
        if kind == 'protein':
            columns.append('in_region')
        columns += ['CAAS_score', 'status']
        summary_columns = ['collection_id', 'term_id', 'n_positions_evaluated',
                           'n_positions_measured', 'coverage_proportion',
                           'median_distance', 'mean_distance', 'n_genes_with_sites']
        summary_columns += [f'prop_within_{r}' for r in radios] + ['status']
        self.detail = self.writer(stack, output_dir, kind, 'distances', columns)
        self.summary = self.writer(stack, output_dir, kind, 'summary', summary_columns)
        self.rows_written = 0
        self.reset()

    @staticmethod
    def writer(stack, output_dir, kind, suffix, columns):
        path = os.path.join(output_dir, f'proximity_{kind}_{suffix}.tsv')
        stream = stack.enter_context(open(path, 'w', newline=''))
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        return writer

    def reset(self):
        self.distances = array('d')
        self.total = 0
        self.distance_sum = 0
        self.within = [0] * len(self.radios)
        self.genes = set()

    def append(self, row):
        self.detail.writerow(row)
        self.rows_written += 1
        self.total += 1
        self.genes.add(row['Gene'])
        if row['status'] == 'measured':
            distance = row[self.dist_col]
            self.distances.append(distance)
            self.distance_sum += distance
            for index, radius in enumerate(self.radios):
                self.within[index] += distance <= radius

    def finish_term(self, cid, term):
        if self.total:
            measured = len(self.distances)
            # Partition the compact buffer in place: no detail dataframe or sort copy.
            median = float(np.median(np.frombuffer(self.distances, dtype=np.float64),
                                     overwrite_input=True)) if measured else None
            row = dict(collection_id=cid, term_id=term, n_positions_evaluated=self.total,
                       n_positions_measured=measured,
                       coverage_proportion=repr(measured / self.total),
                       median_distance=repr(median) if measured else 'NA',
                       mean_distance=repr(self.distance_sum / measured) if measured else 'NA',
                       n_genes_with_sites=len(self.genes),
                       status='descriptive')
            for radius, count in zip(self.radios, self.within):
                row[f'prop_within_{radius}'] = repr(count / measured) if measured else repr('NA')
            self.summary.writerow(row)
        self.reset()


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

    collections = list(read_tsv(os.path.join(args.registry_dir, 'collections.tsv')))
    with open(os.path.join(args.registry_dir, 'manifest.json')) as f:
        manifest = json.load(f)
    wanted = set(args.collections.split(',')) if args.collections else set(PROTEIN_COLLECTIONS)
    selected_collections = {r['collection_id'] for r in collections
                            if r['status'] == 'included' and r['collection_id'] in wanted
                            and r['collection_id'] in PROTEIN_COLLECTIONS}
    region_type = {(r['collection_id'], r['term_id']): r['term_type'] for r in collections}

    # Stream registry inputs instead of materializing dataframes and record copies.
    term_members = {}
    needed = {}
    for row in pm.member_records(args.registry_dir, selected_collections):
        key = (row['collection_id'], row['term_id'])
        gene, pos = row['Gene'], int(row['Position'])
        term_members.setdefault(key, {}).setdefault(gene, set()).add(pos)
        needed.setdefault(gene, set()).add(pos)

    # One best-side score per position, with queries indexed by gene.
    scores = {gene: {} for gene in needed}
    for row in read_tsv(os.path.join(args.registry_dir, 'position_evidence.tsv')):
        gene = row['Gene']
        if gene not in scores or row['CAAS_score'] in ('NA', ''):
            continue
        pos, score = int(row['Position']), float(row['CAAS_score'])
        current = scores[gene].get(pos)
        if current is None or score > current:
            scores[gene][pos] = score
        needed[gene].add(pos)

    coords = {}
    for row in read_tsv(os.path.join(args.registry_dir, 'positions.tsv')):
        gene, pos = row['Gene'], int(row['Position'])
        if gene in needed and pos in needed[gene]:
            coords[(gene, pos)] = {
                'residue': None if row['hg38_aa_pos'] in ('NA', '') else int(row['hg38_aa_pos']),
                'nt': parse_coord(row['hg38_nt_coord']),
            }
    del needed

    with ExitStack() as stack:
        protein_rows = DistanceOutput(stack, args.output_dir, 'protein', 'hg38_aa_pos',
                                      'distance_residues', radios)
        genomic_rows = DistanceOutput(stack, args.output_dir, 'genomic', 'hg38_nt_coord',
                                      'distance_nt', radios)
        near_radius = max(radios, default=0)
        hits = DistanceOutput.writer(stack, args.output_dir, 'protein', 'hits', HIT_COLUMNS)
        by_gene_out = DistanceOutput.writer(stack, args.output_dir, 'protein', 'by_gene', BY_GENE_COLUMNS)
        for (cid, term), by_gene in sorted(term_members.items()):
            is_region = region_type.get((cid, term)) == 'region'
            for gene in sorted(by_gene):
                site_positions = by_gene[gene]
                site_residues = [coords[(gene, p)]['residue'] for p in sorted(site_positions)
                                 if (gene, p) in coords]
                site_residues = sorted(r for r in site_residues if r is not None)
                site_nt = [coords[(gene, p)]['nt'] for p in sorted(site_positions)
                           if (gene, p) in coords]
                site_nt = [c for c in site_nt if c is not None]
                sites_by_chrom = {}
                for chrom, nt in site_nt:
                    sites_by_chrom.setdefault(chrom, []).append(nt)
                sites_by_chrom = {chrom: sorted(values) for chrom, values in sites_by_chrom.items()}
                residue_position = {coords[(gene, p)]['residue']: p for p in sorted(site_positions)
                                    if (gene, p) in coords}
                n_detections = n_overlap = n_near = 0
                closest = None
                for pos, score in sorted(scores[gene].items()):
                    g = gene
                    info = coords.get((g, pos))
                    if info is None:
                        continue
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
                        identical = pos in site_positions
                        if identical:
                            distance, site_residue = 0, info['residue']
                        else:
                            distance, site_residue = nearest_site(info['residue'], site_residues)
                        n_detections += 1
                        closest = distance if closest is None else min(closest, distance)
                        if identical or distance <= near_radius:
                            n_overlap += identical
                            n_near += not identical
                            hits.writerow({'collection_id': cid, 'term_id': term, 'Gene': g, 'Position': pos,
                                           'hg38_aa_pos': info['residue'],
                                           'relation': 'overlap' if identical else 'near',
                                           'distance_residues': distance,
                                           'nearest_site_position': pos if identical else residue_position.get(site_residue, 'NA'),
                                           'nearest_site_residue': site_residue, 'in_region': in_region,
                                           'CAAS_score': score})
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
                        same_chrom = sites_by_chrom.get(info['nt'][0], [])
                        if not same_chrom:
                            genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                                 'Position': pos,
                                                 'hg38_nt_coord': f'{info["nt"][0]}:{info["nt"][1]}',
                                                 'distance_nt': 'NA', 'CAAS_score': score,
                                                 'status': 'no_sites_same_chromosome'})
                        else:
                            distance = 0 if pos in site_positions else nearest_distance(info['nt'][1], same_chrom)
                            genomic_rows.append({'collection_id': cid, 'term_id': term, 'Gene': g,
                                                 'Position': pos,
                                                 'hg38_nt_coord': f'{info["nt"][0]}:{info["nt"][1]}',
                                                 'distance_nt': distance,
                                                 'CAAS_score': score, 'status': 'measured'})

                by_gene_out.writerow({'collection_id': cid, 'term_id': term, 'Gene': gene,
                                      'n_term_positions_in_gene': len(site_positions),
                                      'n_detections_measured': n_detections, 'n_overlaps': n_overlap,
                                      'n_near': n_near,
                                      'min_distance_residues': 'NA' if closest is None else closest})

            protein_rows.finish_term(cid, term)
            genomic_rows.finish_term(cid, term)

    out_manifest = {
        'component': 'proximity',
        'ranking_population': 'canonical_positive_detections_only',
        'contract_version': manifest.get('contract_version'),
        'parameters': {'collections': args.collections, 'radios': radios},
        'contracts': {'protein': 'distancia en residuos de referencia (hg38_aa_pos) dentro de la misma proteina; 0 dentro de region',
                      'genomic': 'distancia en nt dentro del mismo gen y cromosoma; 0 en coincidencia'},
        'states': {'proximity': 'executed_descriptive' if protein_rows.rows_written else 'no_terms',
                   'characterization': 'not_executed (componente separado)'},
        'limits': ['sin p: las distancias son descriptivas',
                   'radios descriptivos fijados antes del contraste',
                   'overlap es identidad exacta de la posicion; near es una posicion distinta dentro del radio mayor',
                   'proteinas sin referencia comparable se informan como falta de cobertura, no distancia infinita',
                   'proximidad genómica y proteica no se mezclan',
                   'la distancia no se mide en columnas MSA con gaps',
                   'consulta condicionada a detecciones positivas; no es el universo completo evaluado'],
    }
    with open(os.path.join(args.output_dir, 'proximity_manifest.json'), 'w') as f:
        json.dump(out_manifest, f, indent=2, sort_keys=True)
        f.write('\n')
    print(f"proximity: {protein_rows.rows_written} protein distance rows, {genomic_rows.rows_written} genomic rows "
          f"-> {args.output_dir}")


if __name__ == '__main__':
    main()
