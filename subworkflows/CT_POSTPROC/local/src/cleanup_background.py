#!/usr/bin/env python3
# cleanup_background.py — Remove the genes dropped by the gene filter from the background gene list.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Background cleanup for CAAS post-processing, from one global background.

1) Load the global background genes (one gene per line, or a TSV with a Gene column).
2) Load removed_genes_summary.tsv from the gene filter.
3) Write cleaned_background_main.txt: the background minus the removed genes.

Called by:  CAAS_BACKGROUND_CLEANUP process (ctpp_clustfilter.nf)
Inputs:     -s  removed_genes_summary.tsv (Gene, category)
            -g  global background gene list
            -o  output directory
Outputs:    cleaned_background_main.txt
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def load_removed_genes(summary_file):
    """Read the removed-genes summary and return its genes as a set."""
    if not Path(summary_file).exists():
        print(f"Error: Removed genes summary file not found: {summary_file}", file=sys.stderr)
        sys.exit(1)

    removed_df = pd.read_csv(summary_file, sep='\t', dtype=str, keep_default_na=False)
    if 'Gene' not in removed_df.columns:
        print("Error: Removed genes summary missing columns: ['Gene']", file=sys.stderr)
        sys.exit(1)

    removed = set(removed_df['Gene'].dropna().astype(str))
    print(f"Loaded {len(removed)} removed genes", file=sys.stderr)
    return removed


def load_global_background_genes(bg_file):
    """Read the global background gene universe: a Gene column when present, else the first column."""
    p = Path(bg_file)
    if not p.exists():
        print(f"Error: Global background file not found: {bg_file}", file=sys.stderr)
        sys.exit(1)

    # A single-column list may be empty and identifiers such as NA are literal.
    import csv
    with p.open(newline='') as fh:
        rows = [row for row in csv.reader(fh, delimiter='\t') if row and not row[0].startswith('#')]
    if rows and 'Gene' in rows[0]:
        col = rows.pop(0).index('Gene')
    else:
        col = 0
    if any(len(row) <= col for row in rows):
        raise ValueError('Malformed background row')
    genes = [row[col] for row in rows]
    genes = sorted({g.strip() for g in genes if g.strip()})
    print(f"Loaded global background genes: {len(genes)}", file=sys.stderr)
    return genes


def write_gene_list(genes, output_file):
    """One gene per line."""
    Path(output_file).write_text("\n".join(genes) + ("\n" if genes else ""))


def main():
    parser = argparse.ArgumentParser(description="Global background cleanup")
    parser.add_argument('-s', '--summary', required=True, help='removed_genes_summary.tsv')
    parser.add_argument('-g', '--global-background-file', required=True,
                        help='Global background genes file (e.g. background_genes.output)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    removed = load_removed_genes(args.summary)
    bg_genes = load_global_background_genes(args.global_background_file)

    cleaned = sorted(set(bg_genes) - removed)
    write_gene_list(cleaned, output_dir / "cleaned_background_main.txt")
    print(f"Background: removed {len(set(bg_genes) & removed)} genes, retained {len(cleaned)}", file=sys.stderr)

    print("\n=== Background Cleanup Complete ===", file=sys.stderr)


if __name__ == '__main__':
    main()
