#!/usr/bin/env python3
"""
CAAS Gene-Level Filtering Module

Implements dual filtering strategy for CAAS discovery results:
1. Extreme genes: Top percentile by CAAS density (n_CAAS/Length)
2. Dubious genes: IQR outliers (Q3 + k*IQR) with clustered positions

Author: PhyloPhere Pipeline
License: GPL-3.0
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np


def detect_extreme_genes(discovery_df, gene_length_df, percentile=0.99, trait_col="Trait"):
    """
    Identify genes in the top percentile by CAAS density.
    
    CAAP-aware: If discovery_df has CAAP_Group column, thresholds are calculated
    independently per group (genes are evaluated within their group context).
    
    Args:
        discovery_df: DataFrame with columns [Trait, Gene, Position, ...] and optional CAAP_Group
        gene_length_df: DataFrame with columns [Gene, Length]
        percentile: Threshold percentile (default 0.99 for top 1%)
        trait_col: Name of trait column
    
    Returns:
        DataFrame with columns [Trait, Gene, n_CAAS, Length, n_CAAS_per_length, Category] and optional CAAP_Group
    """
    # Check if CAAP mode (CAAP_Group column present)
    has_caap_group = 'CAAP_Group' in discovery_df.columns
    groupby_cols = [trait_col, 'CAAP_Group'] if has_caap_group else [trait_col]
    
    # Group-aware thresholds are always calculated per CAAP_Group when present
    discovery_for_removal = discovery_df.copy()
    if has_caap_group:
        print(f"CAAP mode detected: calculating extreme gene thresholds per group", file=sys.stderr)
    
    # Merge gene lengths
    discovery_with_len = discovery_for_removal.merge(
        gene_length_df[['Gene', 'Length']], 
        on='Gene', 
        how='left'
    )
    
    # Check for missing lengths
    missing_len = discovery_with_len['Length'].isna().sum()
    if missing_len > 0:
        print(f"Warning: {missing_len} CAAS positions lack gene length annotation", file=sys.stderr)
        discovery_with_len = discovery_with_len.dropna(subset=['Length'])
    
    # Count unique CAAS positions per gene per trait (and group if CAAP mode)
    groupby_summary = groupby_cols + ['Gene', 'Position']
    gene_summary = (
        discovery_with_len
        .groupby(groupby_summary)
        .first()  # Keep only one occurrence of each position
        .reset_index()
        .groupby(groupby_cols + ['Gene', 'Length'])
        .size()
        .reset_index(name='n_CAAS')
    )
    
    # Calculate CAAS density (percentage)
    gene_summary['n_CAAS_per_length'] = (gene_summary['n_CAAS'] / gene_summary['Length']) * 100
    
    # Calculate threshold per trait (and per group if CAAP mode)
    thresholds = (
        gene_summary
        .groupby(groupby_cols)['n_CAAS_per_length']
        .quantile(percentile)
        .reset_index(name='threshold')
    )
    
    # Flag extreme genes
    gene_summary = gene_summary.merge(thresholds, on=groupby_cols)
    extreme_genes = gene_summary[gene_summary['n_CAAS_per_length'] > gene_summary['threshold']].copy()
    extreme_genes['Category'] = 'Extreme'
    extreme_genes = extreme_genes.drop(columns=['threshold'])
    
    if has_caap_group:
        group_counts = extreme_genes.groupby('CAAP_Group').size()
        print(f"Detected {len(extreme_genes)} extreme genes (top {100*(1-percentile):.1f}% by density per group):", file=sys.stderr)
        for group, count in group_counts.items():
            print(f"  {group}: {count} genes", file=sys.stderr)
    else:
        print(f"Detected {len(extreme_genes)} extreme genes (top {100*(1-percentile):.1f}% by density)", file=sys.stderr)
    
    return extreme_genes


def detect_dubious_genes(discovery_df, gene_length_df, cluster_file, iqr_multiplier=3.0, trait_col="Trait"):
    """
    Identify IQR outlier genes that contain clustered CAAS positions.
    
    CAAP-aware: If discovery_df has CAAP_Group column, IQR thresholds are calculated
    independently per group (genes are evaluated within their group context).
    
    Args:
        discovery_df: DataFrame with columns [Trait, Gene, Position, ...] and optional CAAP_Group
        gene_length_df: DataFrame with columns [Gene, Length]
        cluster_file: Path to cluster filtering output (*.filtered.*.tsv)
        iqr_multiplier: IQR multiplier for outlier threshold (default 3.0)
        trait_col: Name of trait column
    
    Returns:
        DataFrame with columns [Trait, Gene, n_CAAS, Length, n_CAAS_per_length, Category] and optional CAAP_Group
    """
    # Check if CAAP mode (CAAP_Group column present)
    has_caap_group = 'CAAP_Group' in discovery_df.columns
    groupby_cols = [trait_col, 'CAAP_Group'] if has_caap_group else [trait_col]
    
    # Group-aware thresholds are always calculated per CAAP_Group when present
    discovery_for_removal = discovery_df.copy()
    if has_caap_group:
        print(f"CAAP mode detected: calculating IQR thresholds per group", file=sys.stderr)
    
    # Merge gene lengths
    discovery_with_len = discovery_for_removal.merge(
        gene_length_df[['Gene', 'Length']], 
        on='Gene', 
        how='left'
    )
    
    # Check for missing lengths
    missing_len = discovery_with_len['Length'].isna().sum()
    if missing_len > 0:
        print(f"Warning: {missing_len} CAAS positions lack gene length annotation", file=sys.stderr)
        discovery_with_len = discovery_with_len.dropna(subset=['Length'])
    
    # Count unique CAAS positions per gene per trait (and group if CAAP mode)
    groupby_summary = groupby_cols + ['Gene', 'Position']
    gene_summary = (
        discovery_with_len
        .groupby(groupby_summary)
        .first()  # Keep only one occurrence of each position
        .reset_index()
        .groupby(groupby_cols + ['Gene', 'Length'])
        .size()
        .reset_index(name='n_CAAS')
    )
    
    # Calculate per-trait (and per-group if CAAP) IQR thresholds
    def calc_iqr_threshold(series, multiplier):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        return q3 + multiplier * iqr
    
    thresholds = (
        gene_summary
        .groupby(groupby_cols)['n_CAAS']
        .apply(lambda x: calc_iqr_threshold(x, iqr_multiplier))
        .reset_index(name='iqr_threshold')
    )
    
    # Flag IQR outliers
    gene_summary = gene_summary.merge(thresholds, on=groupby_cols)
    iqr_outliers = gene_summary[gene_summary['n_CAAS'] > gene_summary['iqr_threshold']].copy()
    
    if has_caap_group:
        if len(iqr_outliers) > 0:
            group_counts = iqr_outliers.groupby('CAAP_Group').size()
            print(f"Detected {len(iqr_outliers)} IQR outlier genes ({iqr_multiplier}×IQR threshold per group):", file=sys.stderr)
            for group, count in group_counts.items():
                print(f"  {group}: {count} genes", file=sys.stderr)
        else:
            print(f"Detected 0 IQR outlier genes ({iqr_multiplier}×IQR threshold per group)", file=sys.stderr)
    else:
        print(f"Detected {len(iqr_outliers)} IQR outlier genes ({iqr_multiplier}×IQR threshold)", file=sys.stderr)
    
    # Early exit if no IQR outliers found
    if len(iqr_outliers) == 0:
        empty_cols = groupby_cols + ['Gene', 'n_CAAS', 'Length', 'n_CAAS_per_length', 'Category']
        print(f"No dubious genes detected (no IQR outliers to check for clusters)", file=sys.stderr)
        return pd.DataFrame(columns=empty_cols)
    
    # Load cluster filtering results
    if not Path(cluster_file).exists():
        print(f"Warning: Cluster file {cluster_file} not found. Cannot identify dubious genes.", file=sys.stderr)
        empty_cols = groupby_cols + ['Gene', 'n_CAAS', 'Length', 'n_CAAS_per_length', 'Category']
        return pd.DataFrame(columns=empty_cols)
    
    cluster_df = pd.read_csv(cluster_file, sep='\t')
    
    # Check required columns
    required_cols = ['Gene', 'Position', 'ClusteringFlag']
    if not all(col in cluster_df.columns for col in required_cols):
        print(f"Error: Cluster file missing required columns: {required_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Identify genes with discarded (clustered) positions
    # If CAAP mode, only consider clusters within the same CAAP group
    if has_caap_group and 'CAAP_Group' in cluster_df.columns:
        # Build set of (Gene, CAAP_Group) tuples with clusters
        genes_with_clusters = set(
            cluster_df[cluster_df['ClusteringFlag'] == 'Discarded']
            [['Gene', 'CAAP_Group']]
            .itertuples(index=False, name=None)
        )
        # Filter iqr_outliers: keep only those with matching (Gene, CAAP_Group) in clusters
        dubious_genes = iqr_outliers[
            iqr_outliers.apply(lambda row: (row['Gene'], row['CAAP_Group']) in genes_with_clusters, axis=1)
        ].copy()
    else:
        # Standard mode: just check gene presence
        genes_with_clusters = cluster_df[cluster_df['ClusteringFlag'] == 'Discarded']['Gene'].unique()
        dubious_genes = iqr_outliers[iqr_outliers['Gene'].isin(genes_with_clusters)].copy()
    
    # Only modify non-empty DataFrame
    if len(dubious_genes) > 0:
        dubious_genes['Category'] = 'Dubious'
        # Drop iqr_threshold column if it exists
        if 'iqr_threshold' in dubious_genes.columns:
            dubious_genes = dubious_genes.drop(columns=['iqr_threshold'])
        # Calculate CAAS density for output
        dubious_genes['n_CAAS_per_length'] = (dubious_genes['n_CAAS'] / dubious_genes['Length']) * 100
    
    if has_caap_group:
        group_counts = dubious_genes.groupby('CAAP_Group').size()
        print(f"Detected {len(dubious_genes)} dubious genes (IQR outliers with clusters per group):", file=sys.stderr)
        for group, count in group_counts.items():
            print(f"  {group}: {count} genes", file=sys.stderr)
    else:
        print(f"Detected {len(dubious_genes)} dubious genes (IQR outliers with clusters)", file=sys.stderr)
    
    return dubious_genes


def apply_gene_filter(discovery_df, removed_genes_df, mode='both'):
    """
    Remove flagged genes from discovery dataset.
    
    Args:
        discovery_df: Original CAAS discovery DataFrame
        removed_genes_df: DataFrame of genes to remove with Category column
        mode: Filter mode ('none', 'extreme', 'dubious', 'both')
    
    Returns:
        Filtered discovery DataFrame
    """
    if mode == 'none':
        print("Filter mode 'none': No gene-level filtering applied", file=sys.stderr)
        return discovery_df
    
    # Filter removed_genes_df by mode
    if mode == 'extreme':
        genes_to_remove = removed_genes_df[removed_genes_df['Category'] == 'Extreme']
    elif mode == 'dubious':
        genes_to_remove = removed_genes_df[removed_genes_df['Category'] == 'Dubious']
    elif mode == 'both':
        genes_to_remove = removed_genes_df
    else:
        print(f"Error: Invalid filter mode '{mode}'. Use: none, extreme, dubious, both", file=sys.stderr)
        sys.exit(1)
    
    # Remove genes (group-aware when CAAP_Group is present)
    n_before = len(discovery_df)

    if 'CAAP_Group' in discovery_df.columns and 'CAAP_Group' in genes_to_remove.columns:
        remove_keys = set(
            genes_to_remove[['Gene', 'CAAP_Group']]
            .dropna()
            .itertuples(index=False, name=None)
        )
        filtered_df = discovery_df[
            ~discovery_df.apply(lambda row: (row['Gene'], row['CAAP_Group']) in remove_keys, axis=1)
        ].copy()
        n_gene_units = len(remove_keys)
    else:
        gene_list = genes_to_remove['Gene'].unique()
        filtered_df = discovery_df[~discovery_df['Gene'].isin(gene_list)].copy()
        n_gene_units = len(gene_list)

    n_after = len(filtered_df)
    n_removed = n_before - n_after
    
    print(f"Removed {n_removed} CAAS positions from {n_gene_units} gene/group units (mode: {mode})", file=sys.stderr)
    pct = (100 * n_after / n_before) if n_before > 0 else 100.0
    print(f"Remaining: {n_after}/{n_before} positions ({pct:.1f}%)", file=sys.stderr)
    
    return filtered_df


def build_full_gene_stats(discovery_df, gene_length_df, removed_genes_df, extreme_percentile=0.99, iqr_multiplier=3.0, trait_col="Trait"):
    """
    Build complete gene statistics for all genes with thresholds and categories.
    
    Args:
        discovery_df: Original CAAS discovery DataFrame
        gene_length_df: DataFrame with columns [Gene, Length]
        removed_genes_df: DataFrame of flagged genes with Category column
        extreme_percentile: Threshold percentile for extreme genes
        iqr_multiplier: IQR multiplier for dubious gene threshold
        trait_col: Name of trait column
    
    Returns:
        DataFrame with all genes, their stats, thresholds, and categories
    """
    # Check if CAAP mode (CAAP_Group column present)
    has_caap_group = 'CAAP_Group' in discovery_df.columns
    groupby_cols = [trait_col, 'CAAP_Group'] if has_caap_group else [trait_col]
    
    # Merge gene lengths
    discovery_with_len = discovery_df.merge(
        gene_length_df[['Gene', 'Length']], 
        on='Gene', 
        how='left'
    )
    
    # Drop rows without gene length annotation
    discovery_with_len = discovery_with_len.dropna(subset=['Length'])
    
    # Count unique CAAS positions per gene per trait (and group if CAAP mode)
    groupby_summary = groupby_cols + ['Gene', 'Position']
    gene_summary = (
        discovery_with_len
        .groupby(groupby_summary)
        .first()  # Keep only one occurrence of each position
        .reset_index()
        .groupby(groupby_cols + ['Gene', 'Length'])
        .size()
        .reset_index(name='n_CAAS')
    )
    
    # Calculate CAAS density (percentage)
    gene_summary['n_CAAS_per_length'] = (gene_summary['n_CAAS'] / gene_summary['Length']) * 100
    
    # Calculate extreme threshold per trait (and per group if CAAP mode)
    extreme_thresholds = (
        gene_summary
        .groupby(groupby_cols)['n_CAAS_per_length']
        .quantile(extreme_percentile)
        .reset_index(name='threshold_extreme')
    )
    
    # Calculate IQR thresholds per trait (and per group if CAAP mode)
    def calc_iqr_threshold(series, multiplier):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        return q3 + multiplier * iqr
    
    iqr_thresholds = (
        gene_summary
        .groupby(groupby_cols)['n_CAAS']
        .apply(lambda x: calc_iqr_threshold(x, iqr_multiplier))
        .reset_index(name='threshold_dubious')
    )
    
    # Join thresholds to gene summary
    gene_summary = gene_summary.merge(extreme_thresholds, on=groupby_cols)
    gene_summary = gene_summary.merge(iqr_thresholds, on=groupby_cols)
    
    # Initialize all genes as Normal
    gene_summary['Category'] = 'Normal'
    
    # Apply categories from removed_genes_df if available
    if not removed_genes_df.empty:
        # Determine merge columns
        merge_cols = [trait_col, 'CAAP_Group', 'Gene'] if has_caap_group else [trait_col, 'Gene']
        # Only merge if removed_genes_df has the required columns
        if all(col in removed_genes_df.columns for col in merge_cols):
            category_mapping = removed_genes_df[merge_cols + ['Category']].drop_duplicates()
            gene_summary = gene_summary.merge(
                category_mapping, 
                on=merge_cols, 
                how='left', 
                suffixes=('', '_removed')
            )
            # Update categories for flagged genes
            gene_summary['Category'] = gene_summary['Category_removed'].fillna(gene_summary['Category'])
            gene_summary = gene_summary.drop(columns=['Category_removed'])
    
    return gene_summary


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
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -c clusters.tsv \
    --extreme-percentile 0.95 --iqr-multiplier 4.0 -o filtered.tsv
"""
    )
    
    # Input files
    parser.add_argument('-i', '--disambiguation-input', required=True,
                        help='CAAS disambiguation file (TSV with Trait, Gene, Position columns)')
    parser.add_argument('-l', '--gene-ensembl-file', required=True,
                        help='Gene annotation file (TSV with Gene, Chr, Start, End, Strand, Length columns)')
    parser.add_argument('-c', '--cluster-file', default=None,
                        help='Cluster filtering output (required for dubious gene detection)')
    
    # Filter parameters
    parser.add_argument('-m', '--filter-mode', 
                        choices=['none', 'extreme', 'dubious', 'both'],
                        default='both',
                        help='Filtering mode (default: both)')
    parser.add_argument('--extreme-percentile', type=float, default=0.99,
                        help='Density percentile for extreme genes (default: 0.99 = top 1%%)')
    parser.add_argument('--iqr-multiplier', type=float, default=3.0,
                        help='IQR multiplier for dubious gene threshold (default: 3.0)')
    parser.add_argument('--trait-col', default='Trait',
                        help='Name of trait column (default: Trait)')
    
    # Output files
    parser.add_argument('-o', '--output', required=True,
                        help='Filtered discovery output file')
    parser.add_argument('-s', '--summary', default=None,
                        help='Removed genes summary file (optional)')
    parser.add_argument('-g', '--gene-stats-output', default=None,
                        help='Full gene statistics output file (optional)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.disambiguation_input).exists():
        print(f"Error: Disambiguation file not found: {args.disambiguation_input}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.gene_ensembl_file).exists():
        print(f"Error: Gene ensembl file not found: {args.gene_ensembl_file}", file=sys.stderr)
        sys.exit(1)
    
    if args.filter_mode in ['dubious', 'both'] and args.cluster_file is None:
        print("Error: --cluster-file required for 'dubious' or 'both' filter modes", file=sys.stderr)
        sys.exit(1)
    
    # Load data
    print(f"Loading disambiguation data: {args.disambiguation_input}", file=sys.stderr)
    discovery_df = pd.read_csv(args.disambiguation_input, sep='\t')
    
    print(f"Loading gene lengths: {args.gene_ensembl_file}", file=sys.stderr)
    gene_length_df = pd.read_csv(args.gene_ensembl_file, sep='\t')
    
    # Validate required columns - handle both old format and new ensembl format
    # Old format: Gene, Length
    # New format: gene, chr, start, end, strand, length (case insensitive)
    
    # Normalize column names to lowercase for case-insensitive matching
    gene_length_df.columns = gene_length_df.columns.str.lower()
    
    if 'gene' not in gene_length_df.columns:
        print("Error: Gene ensembl file must have 'gene' column", file=sys.stderr)
        sys.exit(1)
    
    if 'length' not in gene_length_df.columns:
        print("Error: Gene ensembl file must have 'length' column", file=sys.stderr)
        sys.exit(1)
    
    # Rename to standard format for processing (capitalize first letter)
    gene_length_df.rename(columns={'gene': 'Gene', 'length': 'Length'}, inplace=True)
    
    print(f"Loaded {len(discovery_df)} CAAS positions across {discovery_df['Gene'].nunique()} genes", file=sys.stderr)
    
    # Detect extreme genes
    removed_genes_list = []
    
    if args.filter_mode in ['extreme', 'both']:
        extreme_genes = detect_extreme_genes(
            discovery_df, 
            gene_length_df, 
            percentile=args.extreme_percentile,
            trait_col=args.trait_col
        )
        removed_genes_list.append(extreme_genes)
    
    # Detect dubious genes
    if args.filter_mode in ['dubious', 'both']:
        dubious_genes = detect_dubious_genes(
            discovery_df,
            gene_length_df,
            args.cluster_file,
            iqr_multiplier=args.iqr_multiplier,
            trait_col=args.trait_col
        )
        removed_genes_list.append(dubious_genes)
    
    # Combine removed genes
    if removed_genes_list:
        removed_genes_df = pd.concat(removed_genes_list, ignore_index=True)
        
        # Determine groupby columns for duplicate detection
        has_caap_group = 'CAAP_Group' in removed_genes_df.columns
        dup_cols = [args.trait_col, 'CAAP_Group', 'Gene'] if has_caap_group else [args.trait_col, 'Gene']
        
        # Handle genes flagged by both criteria (within same group if CAAP mode)
        duplicate_genes = removed_genes_df[removed_genes_df.duplicated(subset=dup_cols, keep=False)]
        if len(duplicate_genes) > 0:
            dup_count = len(duplicate_genes)//2
            if has_caap_group:
                print(f"Found {dup_count} gene-group combinations flagged as both Extreme and Dubious", file=sys.stderr)
            else:
                print(f"Found {dup_count} genes flagged as both Extreme and Dubious", file=sys.stderr)
            # Mark as "Both" and deduplicate
            removed_genes_df.loc[
                removed_genes_df.duplicated(subset=dup_cols, keep=False), 
                'Category'
            ] = 'Both'
            removed_genes_df = removed_genes_df.drop_duplicates(subset=dup_cols, keep='first')
    else:
        # Create empty DataFrame with appropriate columns
        has_caap_group = 'CAAP_Group' in discovery_df.columns
        base_cols = [args.trait_col, 'Gene', 'n_CAAS', 'Length', 'n_CAAS_per_length', 'Category']
        cols = [args.trait_col, 'CAAP_Group'] + base_cols[1:] if has_caap_group else base_cols
        removed_genes_df = pd.DataFrame(columns=cols)
    
    # Apply filtering
    filtered_df = apply_gene_filter(discovery_df, removed_genes_df, mode=args.filter_mode)
    
    # Write outputs
    print(f"Writing filtered discovery: {args.output}", file=sys.stderr)
    filtered_df.to_csv(args.output, sep='\t', index=False)
    
    if args.summary:
        print(f"Writing removed genes summary: {args.summary}", file=sys.stderr)
        removed_genes_df.to_csv(args.summary, sep='\t', index=False)
    
    # Generate and export full gene statistics if requested
    if args.gene_stats_output:
        print(f"Building full gene statistics...", file=sys.stderr)
        gene_stats_df = build_full_gene_stats(
            discovery_df,
            gene_length_df,
            removed_genes_df,
            extreme_percentile=args.extreme_percentile,
            iqr_multiplier=args.iqr_multiplier,
            trait_col=args.trait_col
        )
        print(f"Writing gene statistics: {args.gene_stats_output}", file=sys.stderr)
        gene_stats_df.to_csv(args.gene_stats_output, sep='\t', index=False)
    
    # Summary statistics
    print("\n=== Gene Filtering Summary ===", file=sys.stderr)
    print(f"Filter mode: {args.filter_mode}", file=sys.stderr)
    print(f"Removed genes by category:", file=sys.stderr)
    for category in ['Extreme', 'Dubious', 'Both']:
        count = len(removed_genes_df[removed_genes_df['Category'] == category])
        if count > 0:
            print(f"  {category}: {count}", file=sys.stderr)
    print(f"Total removed genes: {len(removed_genes_df)}", file=sys.stderr)
    pct_str = f"{100*len(filtered_df)/len(discovery_df):.2f}%" if len(discovery_df) > 0 else "N/A"
    print(f"Positions retained: {len(filtered_df)}/{len(discovery_df)} ({pct_str})", file=sys.stderr)


if __name__ == '__main__':
    main()
