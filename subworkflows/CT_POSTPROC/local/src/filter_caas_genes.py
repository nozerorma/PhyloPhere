#!/usr/bin/env python3
"""
CAAS Gene-Level Filtering Module

Filters the pooled CAAS discovery table by gene-level outlier criteria, using the same
implementation (core.postproc) as the permulation null:
1. Extreme genes: density (distinct positions / gene length) above a percentile.
2. Dubious genes: IQR outliers (Q3 + k*IQR) in distinct positions that also carry a
   cluster-train position.

Thresholds are calibrated within each caap_group over the pooled rows (the observed
labeling), with no per-hypothesis grain. Cluster positions are dropped only with
--remove-clusters.

Author: PhyloPhere Pipeline
License: GPL-3.0
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# core.postproc lives with the disambiguation core; it is the single implementation of
# trains and gene removal for the observed and null chains.
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION" / "local"))
from src.core.postproc import GeneUnit, gene_removal, gene_unit_stats, load_gene_lengths  # noqa: E402

# The CAAS table carries categorical amino-acid columns (caas, amino_encoded,
# derived_residues) whose values can legitimately equal NA-sentinel strings --
# "N/A" is Asn on the changed side against Ala -- and pandas' default NA parsing
# would silently blank them (this is exactly how derived_residues was being lost).
# Only a truly empty cell is missing data in these files.
_CAAS_READ_KW = dict(keep_default_na=False, na_values=["", "nan", "NaN"])

# Label of the observed slice in core.postproc units, and the hyp_id shown in gene_stats.tsv.
OBSERVED = "b_0"
POOL_ID = "ALL"


def read_discarded(cluster_file):
    """Cluster file to the set of (Gene, Position, caap_group) flagged Discarded."""
    cluster_df = pd.read_csv(cluster_file, sep="\t", **_CAAS_READ_KW)
    required = ["Gene", "Position", "clustering_flag"]
    if not all(c in cluster_df.columns for c in required):
        print(f"Error: Cluster file missing required columns: {required}", file=sys.stderr)
        sys.exit(1)
    if "caap_group" not in cluster_df.columns:
        cluster_df["caap_group"] = "US"
    d = cluster_df[cluster_df["clustering_flag"] == "Discarded"]
    return set(d[["Gene", "Position", "caap_group"]].itertuples(index=False, name=None))


def build_units(discovery_df, discarded):
    """One core.postproc unit per (caap_group, Gene) of the pooled table."""
    clustered = {(grp, gene) for gene, _pos, grp in discarded}
    counts = discovery_df.groupby(["caap_group", "Gene"])["Position"].nunique()
    return [GeneUnit(OBSERVED, str(grp), str(gene), int(n), (str(grp), str(gene)) in clustered)
            for (grp, gene), n in counts.items()]


def gene_stats_table(stats, removal):
    """gene_stats.tsv: one row per (caap_group, Gene) with a length; category is the removal verdict."""
    rows = []
    for s in stats:
        if s["density"] is None:
            continue
        cat = removal.get((s["labeling"], s["caap_group"], s["Gene"]), "Normal")
        rows.append({
            "hyp_id": POOL_ID, "caap_group": s["caap_group"], "Gene": s["Gene"],
            "length": s["length"], "n_CAAS": s["n_caas"], "n_CAAS_per_length": s["density"],
            "threshold_extreme": s["threshold_extreme"], "threshold_dubious": s["threshold_dubious"],
            "category": cat,
        })
    cols = ["hyp_id", "caap_group", "Gene", "length", "n_CAAS", "n_CAAS_per_length",
            "threshold_extreme", "threshold_dubious", "category"]
    return pd.DataFrame(rows, columns=cols)


def apply_gene_filter(discovery_df, removed_genes_df):
    """Remove the flagged (Gene, caap_group) units from the discovery table."""
    n_before = len(discovery_df)
    remove_keys = set(removed_genes_df[["Gene", "caap_group"]].itertuples(index=False, name=None))
    keys = pd.Series(list(zip(discovery_df["Gene"], discovery_df["caap_group"])), index=discovery_df.index)
    filtered_df = discovery_df[~keys.isin(remove_keys)].copy()
    n_after = len(filtered_df)
    print(f"Removed {n_before - n_after} CAAS positions from {len(remove_keys)} gene/group units", file=sys.stderr)
    pct = (100 * n_after / n_before) if n_before > 0 else 100.0
    print(f"Remaining: {n_after}/{n_before} positions ({pct:.1f}%)", file=sys.stderr)
    return filtered_df


def apply_cluster_filter(discovery_df, discarded):
    """Remove individual CAAS positions flagged as Discarded by cluster filtering."""
    if not discarded:
        print("Cluster filtering: 0 positions flagged as Discarded", file=sys.stderr)
        return discovery_df
    n_before = len(discovery_df)
    keys = pd.Series(list(zip(discovery_df["Gene"], discovery_df["Position"], discovery_df["caap_group"])),
                     index=discovery_df.index)
    filtered_df = discovery_df[~keys.isin(discarded)].copy()
    n_after = len(filtered_df)
    print(f"Cluster filtering: Removed {n_before - n_after} clustered positions "
          f"({len(discarded)} unique cluster events)", file=sys.stderr)
    pct = (100 * n_after / n_before) if n_before > 0 else 100.0
    print(f"Remaining after cluster filter: {n_after}/{n_before} positions ({pct:.1f}%)", file=sys.stderr)
    return filtered_df


def main():
    parser = argparse.ArgumentParser(
        description="Filter CAAS discovery results by gene-level criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove both extreme and dubious genes
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -c clusters.tsv -m both -o filtered.tsv

  # Remove only extreme genes (top 1% density)
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -m extreme -o filtered.tsv

  # Custom thresholds
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -c clusters.tsv \\
    --extreme-percentile 0.95 --iqr-multiplier 4.0 -o filtered.tsv
"""
    )

    # Input files
    parser.add_argument('-i', '--disambiguation-input', required=True,
                        help='Pooled CAAS discovery file (TSV with Gene, Position and, optionally, caap_group)')
    parser.add_argument('-l', '--gene-ensembl-file', required=True,
                        help='Gene annotation file (TSV with gene and length columns)')
    parser.add_argument('-c', '--cluster-file', default=None,
                        help='Cluster filtering output (required for dubious gene detection)')

    # Filter parameters
    parser.add_argument('-m', '--filter-mode',
                        choices=['none', 'extreme', 'dubious', 'both'],
                        default='both',
                        help='Filtering mode (default: both)')
    parser.add_argument('--remove-clusters', action='store_true', default=False,
                        help='Remove positions flagged as Discarded by cluster filtering, decoupled from gene filter mode')
    parser.add_argument('--extreme-percentile', type=float, default=0.99,
                        help='Density percentile for extreme genes (default: 0.99 = top 1%%)')
    parser.add_argument('--iqr-multiplier', type=float, default=3.0,
                        help='IQR multiplier for dubious gene threshold (default: 3.0)')

    # Output files
    parser.add_argument('-o', '--output', required=True,
                        help='Filtered discovery output file')
    parser.add_argument('-s', '--summary', default=None,
                        help='Removed genes summary file (optional)')
    parser.add_argument('-g', '--gene-stats-output', default=None,
                        help='Full gene statistics output file (optional)')

    args = parser.parse_args()

    for label, path in (("Disambiguation file", args.disambiguation_input),
                        ("Gene ensembl file", args.gene_ensembl_file)):
        if not Path(path).exists():
            print(f"Error: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    needs_clusters = args.remove_clusters or args.filter_mode in ['dubious', 'both']
    if needs_clusters and args.cluster_file is None:
        print("Error: --cluster-file required when --remove-clusters is set or for 'dubious'/'both' filter modes", file=sys.stderr)
        sys.exit(1)
    if needs_clusters and not Path(args.cluster_file).exists():
        print(f"Error: Cluster file not found: {args.cluster_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading disambiguation data: {args.disambiguation_input}", file=sys.stderr)
    discovery_df = pd.read_csv(args.disambiguation_input, sep='\t', **_CAAS_READ_KW)
    if 'caap_group' not in discovery_df.columns:
        discovery_df['caap_group'] = 'US'

    gene_lengths = load_gene_lengths(args.gene_ensembl_file)
    if not gene_lengths:
        print("Error: Gene ensembl file must have 'gene' and 'length' columns", file=sys.stderr)
        sys.exit(1)

    discarded = read_discarded(args.cluster_file) if needs_clusters else set()
    print(f"Loaded {len(discovery_df)} CAAS positions across {discovery_df['Gene'].nunique()} genes", file=sys.stderr)

    units = build_units(discovery_df, discarded)
    removal = gene_removal(units, gene_lengths, args.filter_mode, args.iqr_multiplier, args.extreme_percentile)
    removed_genes_df = pd.DataFrame(
        sorted((grp, gene, cat) for (_lab, grp, gene), cat in removal.items()),
        columns=['caap_group', 'Gene', 'category'])

    filtered_df = apply_gene_filter(discovery_df, removed_genes_df)

    if args.remove_clusters:
        print("Applying cluster position removal (--remove-clusters enabled)...", file=sys.stderr)
        filtered_df = apply_cluster_filter(filtered_df, discarded)
    else:
        print("Cluster position removal: disabled (cluster positions retained unless removed by gene filter)", file=sys.stderr)

    print(f"Writing filtered discovery: {args.output}", file=sys.stderr)
    filtered_df.to_csv(args.output, sep='\t', index=False)

    if args.summary:
        print(f"Writing removed genes summary: {args.summary}", file=sys.stderr)
        removed_genes_df.to_csv(args.summary, sep='\t', index=False)

    if args.gene_stats_output:
        stats = gene_unit_stats(units, gene_lengths, args.iqr_multiplier, args.extreme_percentile)
        print(f"Writing gene statistics: {args.gene_stats_output}", file=sys.stderr)
        gene_stats_table(stats, removal).to_csv(args.gene_stats_output, sep='\t', index=False)

    print("\n=== Post-Processing Filtering Summary ===", file=sys.stderr)
    print(f"Gene filter mode: {args.filter_mode}", file=sys.stderr)
    print(f"Cluster removal: {'enabled' if args.remove_clusters else 'disabled'}", file=sys.stderr)
    for category in ['Extreme', 'Dubious', 'Both']:
        count = int((removed_genes_df['category'] == category).sum())
        if count > 0:
            print(f"  {category}: {count}", file=sys.stderr)
    print(f"Total removed gene/scheme units: {len(removed_genes_df)}", file=sys.stderr)
    pct_str = f"{100*len(filtered_df)/len(discovery_df):.2f}%" if len(discovery_df) > 0 else "N/A"
    print(f"Positions retained: {len(filtered_df)}/{len(discovery_df)} ({pct_str})", file=sys.stderr)


if __name__ == '__main__':
    main()
