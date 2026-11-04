#!/usr/bin/env python3
"""Project scored alignment columns onto the declared MAP nucleotide anchors."""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scores', required=True)
    parser.add_argument('--map-dir', required=True)
    parser.add_argument('--vep-src', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    sys.path.insert(0, args.vep_src)
    import vep_common as vr

    directory = Path(args.map_dir)
    if not directory.is_dir():
        raise ValueError('MAP input must be a directory')
    vr.require(args.scores, ['Gene', 'Position'])
    keys = set()
    for row in vr.table(args.scores):
        gene = row['Gene']
        if not gene or gene != gene.strip() or '/' in gene or '\\' in gene or gene in {'.', '..'}:
            raise ValueError(f'Invalid literal gene ID {gene!r}')
        keys.add((gene, vr.integer(row['Position'], 0)))
    paths = {p.name.removesuffix('.map.tsv'): p for p in directory.glob('*.map.tsv') if p.is_file()}
    genes = sorted({gene for gene, _ in keys})
    resolved = {gene: vr.resolve_gene(gene, paths) for gene in genes}
    used = {}
    for gene, match in resolved.items():
        if match in paths:
            if match in used and used[match] != gene:
                raise ValueError(f'Ambiguous producer gene alias: {match} matches {used[match]}, {gene}')
            used[match] = gene
    maps = {gene: vr.load_map(paths[match]) for gene, match in resolved.items() if match in paths}
    def rows():
        for gene, position in sorted(keys):
            mapping = maps.get(gene)
            entry = mapping.get(position) if mapping is not None else None
            coordinate = entry.get('anchor') if entry else None
            status = ('no_map' if mapping is None else 'position_not_in_map' if entry is None else
                      'no_genomic_coordinate' if coordinate is None else 'mapped')
            yield dict(Gene=gene, Position=position, hg38_nt_coord=coordinate or '',
                       mapping_status=status, map_gene=resolved[gene] if mapping is not None else '')
    vr.write_table(args.output, ['Gene', 'Position', 'hg38_nt_coord', 'mapping_status', 'map_gene'], rows())


if __name__ == '__main__':
    main()
