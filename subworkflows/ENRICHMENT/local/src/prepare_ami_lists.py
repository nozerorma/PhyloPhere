#!/usr/bin/env python3
"""Canonical DOMINO/report seeds. No inference or module computation is performed."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_genes(path):
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def prepare(method, sources, output):
    files = []
    for source in sources:
        p = Path(source)
        if p.is_dir():
            files.extend(p.glob('active_*.tsv' if method == 'caas' else '*.txt'))
        elif p.is_file() and not p.name.startswith('NO_'):
            files.append(p)
    by_name = {}
    for p in files:
        if p.name in by_name and p.read_bytes() != by_name[p.name].read_bytes():
            raise ValueError(f'Conflicting seed sources for {p.name}')
        by_name[p.name] = p
    lists = {}
    if method == 'caas':
        for rank in ('global', 'top', 'bottom'):
            name = f'active_{rank}.tsv'
            if name not in by_name:
                raise ValueError(f'Missing canonical CAAS seed table: {name}')
            with by_name[name].open() as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                if not {'Gene', 'ranking', 'nominal'}.issubset(reader.fieldnames or []):
                    raise ValueError(f'{name}: canonical Gene/ranking/nominal columns required')
                rows = list(reader)
            if any(row['ranking'] != rank or row['nominal'].lower() not in ('true', '1') for row in rows):
                raise ValueError(f'{name}: incompatible ranking or non-active rows')
            genes = [row['Gene'].strip() for row in rows]
            if any(not g for g in genes) or len(set(genes)) != len(genes):
                raise ValueError(f'{name}: empty or duplicate seed gene')
            lists[rank] = set(genes)
    else:
        wanted = ('fade_top_significant', 'fade_bottom_significant') if method == 'fade' else ('rer_significant', 'rer_accelerating', 'rer_decelerating')
        for name in wanted:
            path = by_name.get(name + '.txt')
            if path is not None:
                lists[name] = read_genes(path)
        if method == 'fade' and lists:
            union = set().union(*lists.values())
            old = by_name.get('fade_global_significant.txt')
            if old is not None and read_genes(old) != union:
                raise ValueError('FADE global seed list differs from the union of available sides')
            lists['fade_global_significant'] = union
        if not lists or (method == 'rer' and 'rer_significant' not in lists):
            raise ValueError(f'Missing canonical {method.upper()} seed lists')
        if method == 'rer':
            positive = lists.get('rer_accelerating', set())
            negative = lists.get('rer_decelerating', set())
            if positive & negative:
                raise ValueError('RER positive and negative seed lists overlap')
            if not (positive | negative) <= lists['rer_significant']:
                raise ValueError('Directional RER seeds are absent from the global RER list')
    output.mkdir(parents=True, exist_ok=True)
    # Reusing a local regeneration scratch directory must not retain old view lists.
    for p in output.glob('*.txt'):
        p.unlink()
    manifest = {'contract_version': 1, 'method': method, 'sources': [str(p) for p in sorted(files)], 'lists': {}}
    for name, genes in sorted(lists.items()):
        data = ''.join(g + '\n' for g in sorted(genes)).encode()
        (output / (name + '.txt')).write_bytes(data)
        manifest['lists'][name] = {'n_genes': len(genes), 'sha256': hashlib.sha256(data).hexdigest()}
    (output / 'canonical_lists.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', choices=('caas', 'fade', 'rer'), required=True)
    parser.add_argument('--sources', nargs='+', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.method, args.sources, args.output_dir)
