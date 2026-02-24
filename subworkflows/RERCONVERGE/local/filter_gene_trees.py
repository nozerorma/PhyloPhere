#!/usr/bin/env python3
"""
filter_gene_trees.py
────────────────────
Filters a RERConverge gene-trees file to keep only genes present in a
gene-list file.

Gene-trees format expected (tab-separated, one gene per line):
    GENE_NAME\tNewick_tree_string

Gene-list format: one gene name per line; empty lines and # comments ignored.

Usage:
    python filter_gene_trees.py <gene_trees_file> <gene_list_file> <output_file>

Arguments:
    gene_trees_file   Input gene trees file (tab-separated GENE\tNewick).
    gene_list_file    Plain-text gene list (one gene per line).
    output_file       Output filtered gene trees file.

Author: Miguel Ramon (miguel.ramon@upf.edu)
"""

import sys
import os


def read_gene_list(path: str) -> set:
    """Read gene names from a plain-text file (one per line, # comments ignored)."""
    genes = set()
    with open(path) as fh:
        for line in fh:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    return genes


def filter_gene_trees(input_path: str, gene_list: set, output_path: str):
    """
    Stream through the gene-trees file, writing only lines whose gene name
    (first tab-separated field) appears in gene_list.
    """
    kept = 0
    total = 0
    with open(input_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            line_stripped = line.rstrip('\n')
            if not line_stripped:
                continue
            # Split on the first tab; gene name is field 0
            parts = line_stripped.split('\t', 1)
            gene_name = parts[0].strip()
            total += 1
            if gene_name in gene_list:
                fout.write(line_stripped + '\n')
                kept += 1

    print(
        f"filter_gene_trees: kept {kept} / {total} genes "
        f"({len(gene_list)} in gene list)",
        file=sys.stderr
    )

    missing = gene_list - {
        line.split('\t', 1)[0].strip()
        for line in open(input_path)
        if line.strip()
    }
    if missing:
        print(
            f"filter_gene_trees: WARNING — {len(missing)} genes in list "
            f"not found in gene trees file (will be silently absent): "
            f"{', '.join(sorted(missing)[:10])}"
            f"{'...' if len(missing) > 10 else ''}",
            file=sys.stderr
        )


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: filter_gene_trees.py <gene_trees_file> <gene_list_file> <output_file>",
            file=sys.stderr
        )
        sys.exit(1)

    gene_trees_file = sys.argv[1]
    gene_list_file  = sys.argv[2]
    output_file     = sys.argv[3]

    for path in (gene_trees_file, gene_list_file):
        if not os.path.exists(path):
            print(f"Error: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    gene_list = read_gene_list(gene_list_file)
    if not gene_list:
        print("Warning: gene list is empty — output will be empty.", file=sys.stderr)

    filter_gene_trees(gene_trees_file, gene_list, output_file)


if __name__ == '__main__':
    main()
