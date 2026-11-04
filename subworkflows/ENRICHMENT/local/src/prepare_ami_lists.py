#!/usr/bin/env python3
"""Canonical DOMINO/report seeds. No inference or module computation is performed."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

# Level lists of the BH classes: active_<ranking>_q<level>.tsv (CAAS) and rer_<kind>_q<level>.txt (RER).
CAAS_LEVEL = re.compile(r'^active_(global|top|bottom)_q([0-9.]+)\.tsv$')
RER_LEVEL = re.compile(r'^(rer_significant|rer_accelerating|rer_decelerating)_q([0-9.]+)\.txt$')
FAMILY = re.compile(r'^(.*)_q([0-9.]+)$')


def below(label, level):
    """BH class label '<0.05' is below `level` when its own level is not above it; 'ns' and 'NA' never are."""
    return label.startswith('<') and float(label[1:]) <= level


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
                reader = csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
                if not {'Gene', 'ranking', 'nominal'}.issubset(reader.fieldnames or []):
                    raise ValueError(f'{name}: canonical Gene/ranking/nominal columns required')
                rows = list(reader)
            if any(row['ranking'] != rank or row['nominal'].lower() not in ('true', '1') for row in rows):
                raise ValueError(f'{name}: incompatible ranking or non-active rows')
            genes = [row['Gene'] for row in rows]
            if any(g!=g.strip() for g in genes):raise ValueError('Whitespace in canonical gene identifier')
            if any(not g for g in genes) or len(set(genes)) != len(genes):
                raise ValueError(f'{name}: empty or duplicate seed gene')
            lists[rank] = set(genes)
        for name, path in sorted(by_name.items()):
            m = CAAS_LEVEL.match(name)
            if not m:
                continue
            rank, token = m.groups()
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
                if not {'Gene', 'ranking', 'bh_class'}.issubset(reader.fieldnames or []):
                    raise ValueError(f'{name}: canonical Gene/ranking/bh_class columns required')
                rows = list(reader)
            if any(row['ranking'] != rank or not below(row['bh_class'], float(token)) for row in rows):
                raise ValueError(f'{name}: incompatible ranking or genes not below the BH level')
            genes = [row['Gene'] for row in rows]
            if any(not g or g != g.strip() for g in genes) or len(set(genes)) != len(genes):
                raise ValueError(f'{name}: empty, padded or duplicate seed gene')
            lists[f'{rank}_q{token}'] = set(genes)
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
            for name, path in sorted(by_name.items()):
                m = RER_LEVEL.match(name)
                if m:
                    lists[f'{m.group(1)}_q{m.group(2)}'] = read_genes(path)
            for token in sorted({RER_LEVEL.match(n).group(2) for n in by_name if RER_LEVEL.match(n)}):
                tag = f'_q{token}'
                total, pos, neg = (lists.get(k + tag) for k in ('rer_significant', 'rer_accelerating', 'rer_decelerating'))
                if total is None or pos is None or neg is None:
                    raise ValueError(f'RER level q{token}: the global, positive and negative lists are all required')
                if pos & neg or not (pos | neg) <= total:
                    raise ValueError(f'RER level q{token}: directional seeds overlap or are absent from the global list')
    # A level list equal to the previous level of its family has the same DOMINO result, and an empty one has no seeds: neither is
    # written (DOMINO would run twice or on nothing); the manifest records both.
    families = {}
    for name in lists:
        m = FAMILY.match(name)
        if m:
            families.setdefault(m.group(1), []).append((float(m.group(2)), name))
    aliases, empty = {}, []
    for members in families.values():
        previous = None
        for _, name in sorted(members):
            if not lists[name]:
                empty.append(name)
                del lists[name]
            elif previous is not None and lists[name] == lists[previous]:
                aliases[name] = previous
                del lists[name]
            else:
                previous = name
    output.mkdir(parents=True, exist_ok=True)
    # Reusing a local regeneration scratch directory must not retain old view lists.
    for p in output.glob('*.txt'):
        p.unlink()
    manifest = {'contract_version': 1, 'method': method, 'sources': [str(p) for p in sorted(files)], 'lists': {},
                'aliases': aliases, 'empty': sorted(empty)}
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
