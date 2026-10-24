#!/usr/bin/env python3
"""Recover explicit CAAS columns from legacy Pfam hits without rerunning hmmscan.

Requires the original hmmscan query FASTA and domtblout, plus the productive
trimmed CAAS alignments. Each legacy hit must occur in domtblout and each query
must exactly match the reference sequence after removing gaps/X. Incompatible
or ambiguous inputs abort before creating outputs. Original files are untouched.
Run substantial migrations on a compute node, not the cluster login node.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from map_domain_variability import read_fasta


REQUIRED = {'gene', 'pfam_id', 'target_name', 'description', 'clan_acc',
            'clan_name', 'ali_start', 'ali_end'}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def gene_key(value):
    return value.split('.')[0]


def unique_index(items, label):
    result = {}
    for name, value in items:
        key = gene_key(name)
        if key in result:
            raise ValueError(f'ambiguous {label} for gene {key}')
        result[key] = (name, value)
    return result


def migrate(domains, query_fasta, domtblout, alignment_dir, output,
            ref_species='Homo_sapiens'):
    domains, query_fasta, domtblout, alignment_dir, output = map(
        Path, (domains, query_fasta, domtblout, alignment_dir, output))
    manifest_path = Path(str(output) + '.manifest.json')
    if output.exists() or manifest_path.exists():
        raise ValueError('output or manifest already exists; choose a new destination')
    with domains.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        if REQUIRED - set(fields):
            raise ValueError(f'missing Pfam columns: {sorted(REQUIRED - set(fields))}')
        if 'coordinate_space' in fields or 'alignment_positions' in fields:
            raise ValueError('input already declares coordinates; legacy migration refused')
        rows = list(reader)
    genes = {gene_key(row['gene']) for row in rows}
    queries = unique_index(
        ((header.split()[0], seq) for header, seq in read_fasta(query_fasta)
         if gene_key(header.split()[0]) in genes), 'query')
    alignments = unique_index(
        ((p.stem, p) for p in sorted(alignment_dir.iterdir())
         if p.is_file() and p.suffix.lower() in ('.fa', '.fasta', '.faa', '.fas')
         and gene_key(p.stem) in genes), 'alignment')
    hits = set()
    with domtblout.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 22:
                raise ValueError('malformed hmmscan domtblout row')
            hits.add((parts[3], parts[1].split('.')[0], int(parts[17]), int(parts[18])))
    mappings, audit = {}, []
    for gene in sorted(genes):
        if gene not in queries or gene not in alignments:
            raise ValueError(f'missing original query or productive alignment for {gene}')
        query_id, query = queries[gene]
        _, alignment = alignments[gene]
        refs = [(header, seq) for header, seq in read_fasta(alignment)
                if ref_species.lower() in header.lower()]
        if len(refs) != 1:
            raise ValueError(f'expected one reference sequence for {gene}, got {len(refs)}')
        header, aligned = refs[0]
        columns = [i for i, aa in enumerate(aligned) if aa not in ('-', 'X')]
        ungapped = ''.join(aligned[i] for i in columns)
        if query != ungapped:
            raise ValueError(f'original query differs from productive reference for {gene}; regenerate Pfam')
        mappings[gene] = (query_id, columns)
        audit.append({'gene': gene, 'query_id': query_id, 'reference_header': header,
                      'alignment': str(alignment.resolve()), 'alignment_sha256': sha256(alignment),
                      'query_sequence_sha256': hashlib.sha256(query.encode()).hexdigest()})
    for row in rows:
        gene = gene_key(row['gene'])
        query_id, columns = mappings[gene]
        start, end = int(row['ali_start']), int(row['ali_end'])
        if not 1 <= start <= end <= len(columns):
            raise ValueError(f'Pfam endpoints outside verified query for {gene}: {start}-{end}')
        hit = (query_id, row['pfam_id'].split('.')[0], start, end)
        if hit not in hits:
            raise ValueError(f'legacy hit absent from original hmmscan domtblout: {hit}')
        row['alignment_positions'] = ','.join(map(str, columns[start - 1:end]))
        row['coordinate_space'] = 'CAAS_0based'
    manifest = {'operation': 'verified_legacy_pfam_query_to_caas',
                'source_coordinate_space': 'hmmscan_query_1based',
                'output_coordinate_space': 'CAAS_0based', 'n_hits': len(rows),
                'n_genes': len(genes), 'references': audit,
                'inputs': {str(p.resolve()): sha256(p) for p in (domains, query_fasta, domtblout)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', lineterminator='\n',
                                fieldnames=fields + ['alignment_positions', 'coordinate_space'])
        writer.writeheader()
        writer.writerows(rows)
    manifest['output_sha256'] = sha256(output)
    with manifest_path.open('x') as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--domains', required=True)
    parser.add_argument('--query-fasta', required=True, help='original hmmscan query FASTA')
    parser.add_argument('--domtblout', required=True, help='original hmmscan --domtblout file')
    parser.add_argument('--alignment-dir', required=True, help='productive trimmed CAAS FASTA alignments')
    parser.add_argument('--ref-species', default='Homo_sapiens')
    parser.add_argument('--output', required=True, help='new TSV path; never overwrites')
    args = parser.parse_args()
    try:
        result = migrate(args.domains, args.query_fasta, args.domtblout,
                         args.alignment_dir, args.output, args.ref_species)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'Pfam migration refused: {exc}\n')
    print(f"Verified {result['n_hits']} Pfam hits across {result['n_genes']} genes -> {args.output}")


if __name__ == '__main__':
    main()
