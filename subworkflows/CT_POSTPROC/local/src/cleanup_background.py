#!/usr/bin/env python3
# cleanup_background.py — Remove the genes dropped by the gene filter from the background gene list.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Background cleanup for CAAS post-processing, by caap_group, from one global background.

1) Load the global background genes (one gene per line, or a TSV with a Gene column).
2) Load removed_genes_summary.tsv from the gene filter.
3) Write one cleaned background per caap_group: the background minus the genes removed in it.
4) Write cleaned_background_main.txt: the background minus the genes removed in all the
   groups listed in the summary. A gene removed in only some groups stays in the main list.

Called by:  CAAS_BACKGROUND_CLEANUP process (ctpp_clustfilter.nf)
Inputs:     -s  removed_genes_summary.tsv (Gene, caap_group, category)
            -g  global background gene list
            -o  output directory
Outputs:    cleaned_background_<caap_group>.txt per group and cleaned_background_main.txt
"""

import argparse
import re
import sys
from pathlib import Path
import pandas as pd


def load_removed_genes(summary_file):
    """Read the removed-genes summary; a missing caap_group column is filled with "US"."""
    if not Path(summary_file).exists():
        print(f"Error: Removed genes summary file not found: {summary_file}", file=sys.stderr)
        sys.exit(1)

    removed_df = pd.read_csv(summary_file, sep='\t')

    required_cols = ['Gene']
    missing_cols = [col for col in required_cols if col not in removed_df.columns]
    if missing_cols:
        print(f"Error: Removed genes summary missing columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    if 'caap_group' not in removed_df.columns:
        removed_df['caap_group'] = 'US'
    else:
        removed_df['caap_group'] = removed_df['caap_group'].fillna('US').astype(str)

    print(
        f"Loaded {len(removed_df)} removed gene entries across "
        f"{removed_df['caap_group'].nunique()} groups",
        file=sys.stderr,
    )

    return removed_df


def load_global_background_genes(bg_file):
    """Read the global background gene universe: a Gene column when present, else the first column."""
    p = Path(bg_file)
    if not p.exists():
        print(f"Error: Global background file not found: {bg_file}", file=sys.stderr)
        sys.exit(1)

    try:
        bg_df = pd.read_csv(p, sep='\t', comment='#')
        if 'Gene' in bg_df.columns:
            genes = bg_df['Gene'].dropna().astype(str)
        else:
            bg_df = pd.read_csv(p, sep='\t', header=None, names=['Gene'], comment='#')
            genes = bg_df['Gene'].dropna().astype(str)
    except Exception:
        bg_df = pd.read_csv(p, header=None, names=['Gene'], comment='#')
        genes = bg_df['Gene'].dropna().astype(str)

    genes = sorted({g.strip() for g in genes if g.strip()})
    print(f"Loaded global background genes: {len(genes)}", file=sys.stderr)
    return genes


def sanitize_group_name(group):
    """Group name made safe for a file name (characters outside A-Za-z0-9_.- become underscores)."""
    g = str(group).strip()
    g = re.sub(r'[^A-Za-z0-9_.-]+', '_', g)
    return g or 'UNKNOWN'


def write_gene_list(genes, output_file):
    """One gene per line."""
    Path(output_file).write_text("\n".join(genes) + ("\n" if genes else ""))


def main():
    parser = argparse.ArgumentParser(description="CAAP-aware global background cleanup")
    parser.add_argument('-s', '--summary', required=True, help='removed_genes_summary.tsv')
    parser.add_argument('-g', '--global-background-file', required=True,
                        help='Global background genes file (e.g. background_genes.output)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    removed_genes_df = load_removed_genes(args.summary)
    bg_genes = load_global_background_genes(args.global_background_file)
    bg_gene_set = set(bg_genes)

    groups = sorted(set(removed_genes_df['caap_group'].astype(str)))
    if not groups:
        groups = ['US']

    removed_per_group = {}
    print(f"Detected groups: {groups}", file=sys.stderr)

    for g in groups:
        removed_g = set(
            removed_genes_df.loc[removed_genes_df['caap_group'] == g, 'Gene']
            .dropna()
            .astype(str)
        )
        removed_per_group[g] = removed_g

        cleaned_g = sorted(bg_gene_set - removed_g)
        g_safe = sanitize_group_name(g)
        out_g = output_dir / f"cleaned_background_{g_safe}.txt"
        write_gene_list(cleaned_g, out_g)
        print(
            f"Group {g}: removed {len(removed_g)} genes, retained {len(cleaned_g)}",
            file=sys.stderr,
        )

    # Main list: only genes removed in every group are dropped
    intersection_removed = set.intersection(*removed_per_group.values()) if removed_per_group else set()
    cleaned_main = sorted(bg_gene_set - intersection_removed)
    out_main = output_dir / "cleaned_background_main.txt"
    write_gene_list(cleaned_main, out_main)
    print(
        f"Global(main): removed {len(intersection_removed)} genes present in all groups, "
        f"retained {len(cleaned_main)}",
        file=sys.stderr,
    )

    print("\n=== Background Cleanup Complete ===", file=sys.stderr)


if __name__ == '__main__':
    main()
