#!/usr/bin/env python3
# filter_caas_genes.py — Remove outlier genes (length-stratified fences), and optionally clustered positions, from the pooled CAAS table.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Gene-level filter of the pooled CAAS discovery table, with the implementation of the
permulation null (core.postproc):

A gene is an outlier when its count of distinct positions exceeds the fences of its gene-length stratum
(core.postproc.stratified_fences): Extreme above the far fence (k_far), Dubious above the near fence (k_near) with a
cluster-train position, Both when both hold.

Thresholds are calibrated over the pooled rows (one unit per gene). Cluster positions are dropped only with --remove-clusters.

Called by:  CAAS_FILTER_GENES process (ctpp_clustfilter.nf)
Inputs:     -i  pooled discovery TSV (Gene, Position)
            -l  gene annotation TSV with gene and length columns
            -c  cluster file of filter_caas_clusters-param.py (needed for the dubious and
                both modes and for --remove-clusters)
Outputs:    -o  filtered discovery TSV
            -s  removed (Gene, category) units
            -g  statistics of every candidate gene, including missing lengths (gene_stats.tsv)
"""

from postproc_io import read_table, positions, validate_cluster
import math
import argparse
import sys
from pathlib import Path

import pandas as pd

# core.postproc lives with the disambiguation core and is shared with the null chain.
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "CT_DISAMBIGUATION" / "local"))
from src.core.postproc import GeneUnit, gene_removal, gene_unit_stats, load_gene_lengths  # noqa: E402

# Labeling of the observed slice in core.postproc units, and the hyp_id column of gene_stats.tsv.
OBSERVED = "b_0"
POOL_ID = "ALL"


def read_discarded(cluster_file, discovery=None):
    """Cluster file to the set of (Gene, Position) flagged Discarded."""
    cluster_df = read_table(cluster_file)
    if discovery is not None:
        cluster_df = validate_cluster(discovery, cluster_df)
    required = ["Gene", "Position", "clustering_flag"]
    if not all(c in cluster_df.columns for c in required):
        print(f"Error: Cluster file missing required columns: {required}", file=sys.stderr)
        sys.exit(1)
    d = cluster_df[cluster_df["clustering_flag"] == "Discarded"]
    return set(d[["Gene", "Position"]].itertuples(index=False, name=None))


def build_units(discovery_df, discarded):
    """One core.postproc unit per Gene of the pooled table."""
    clustered = {gene for gene, _pos in discarded}
    counts = discovery_df.groupby("Gene")["Position"].nunique()
    return [GeneUnit(OBSERVED, str(gene), int(n), str(gene) in clustered)
            for gene, n in counts.items()]


def gene_stats_table(stats, removal, length_sources=None):
    """gene_stats.tsv: one row per candidate Gene; category is the removal verdict (Extreme, Dubious, Both or Normal)."""
    rows = []
    for s in stats:
        cat = removal.get((s["labeling"], s["Gene"]), "Normal")
        rows.append({
            "hyp_id": POOL_ID, "Gene": s["Gene"],
            "length_source": (length_sources or {}).get(s["Gene"], "missing"),
            "length": s["length"], "n_CAAS": s["n_caas"], "n_CAAS_per_length": s["density"],
            "threshold_extreme": s["threshold_extreme"], "threshold_dubious": s["threshold_dubious"],
            "category": cat,
        })
    cols = ["hyp_id", "Gene", "length_source", "length", "n_CAAS", "n_CAAS_per_length",
            "threshold_extreme", "threshold_dubious", "category"]
    return pd.DataFrame(rows, columns=cols)


def apply_gene_filter(discovery_df, removed_genes_df):
    """Remove the flagged genes from the discovery table."""
    n_before = len(discovery_df)
    remove_genes = set(removed_genes_df["Gene"])
    filtered_df = discovery_df[~discovery_df["Gene"].isin(remove_genes)].copy()
    n_after = len(filtered_df)
    print(f"Removed {n_before - n_after} CAAS positions from {len(remove_genes)} genes", file=sys.stderr)
    pct = (100 * n_after / n_before) if n_before > 0 else 100.0
    print(f"Remaining: {n_after}/{n_before} positions ({pct:.1f}%)", file=sys.stderr)
    return filtered_df


def apply_cluster_filter(discovery_df, discarded):
    """Remove individual CAAS positions flagged as Discarded by cluster filtering."""
    if not discarded:
        print("Cluster filtering: 0 positions flagged as Discarded", file=sys.stderr)
        return discovery_df
    n_before = len(discovery_df)
    keys = pd.Series(list(zip(discovery_df["Gene"], discovery_df["Position"])), index=discovery_df.index)
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
  # Remove outlier genes with the default fences (10 length strata, k_far 3, k_near 1.5)
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -c clusters.tsv -m stratified -o filtered.tsv

  # Custom fences
  python filter_caas_genes.py -i discovery.tsv -l gene_lengths.tsv -c clusters.tsv \\
    --strata 20 --k-far 4.0 --k-near 2.0 -o filtered.tsv
"""
    )

    parser.add_argument('-i', '--disambiguation-input', required=True,
                        help='Pooled CAAS discovery file (TSV with Gene and Position)')
    parser.add_argument('-l', '--gene-ensembl-file', default=None,
                        help='Gene annotation file (TSV with gene and length columns), used as fallback')
    parser.add_argument('--map-dir', type=str, default=None,
                        help="Directory of the trimmer's per-gene MAP tables for untrimmed lengths")
    parser.add_argument('--map-suffix', type=str, default=".map.tsv",
                        help='File-name tail of the MAP tables')
    parser.add_argument('-c', '--cluster-file', default=None,
                        help='Cluster filtering output (required for dubious gene detection)')

    parser.add_argument('-m', '--filter-mode',
                        choices=['none', 'stratified'],
                        default='stratified',
                        help='Filtering mode (default: stratified)')
    parser.add_argument('--remove-clusters', action='store_true', default=False,
                        help='Remove positions flagged as Discarded by cluster filtering, decoupled from gene filter mode')
    parser.add_argument('--strata', type=int, default=10,
                        help='Number of equal-size gene-length strata (default: 10)')
    parser.add_argument('--k-far', type=float, default=3.0,
                        help='IQR multiplier of the far fence, above which a gene is removed (default: 3.0)')
    parser.add_argument('--k-near', type=float, default=1.5,
                        help='IQR multiplier of the near fence, above which a gene with a cluster train is removed (default: 1.5)')

    parser.add_argument('-o', '--output', required=True,
                        help='Filtered discovery output file')
    parser.add_argument('-s', '--summary', default=None,
                        help='Removed genes summary file (optional)')
    parser.add_argument('-g', '--gene-stats-output', default=None,
                        help='Full gene statistics output file (optional)')

    parser.add_argument('--quality-universe')
    parser.add_argument('--cluster-contract')
    args = parser.parse_args()
    if (any(not math.isfinite(x) or x < 0 for x in [args.k_far, args.k_near])
            or args.strata < 1 or args.k_far < args.k_near):
        parser.error('Invalid IQR multipliers or length strata (k-far must not be below k-near)')

    if not Path(args.disambiguation_input).exists():
        print(f"Error: Disambiguation file not found: {args.disambiguation_input}", file=sys.stderr)
        sys.exit(1)
    if args.gene_ensembl_file and not Path(args.gene_ensembl_file).exists():
        print(f"Error: Gene ensembl file not found: {args.gene_ensembl_file}", file=sys.stderr)
        sys.exit(1)

    needs_clusters = args.remove_clusters or args.filter_mode == 'stratified'
    if needs_clusters and args.cluster_file is None:
        print("Error: --cluster-file required when --remove-clusters is set or for the 'stratified' filter mode", file=sys.stderr)
        sys.exit(1)
    if needs_clusters and not Path(args.cluster_file).exists():
        print(f"Error: Cluster file not found: {args.cluster_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading disambiguation data: {args.disambiguation_input}", file=sys.stderr)
    discovery_df = read_table(args.disambiguation_input)
    discovery_df['Position'] = positions(discovery_df.Position)

    length_sources = {}
    gene_lengths = load_gene_lengths(
        path=args.gene_ensembl_file,
        map_dir=args.map_dir,
        map_suffix=args.map_suffix,
        genes=discovery_df['Gene'].unique(),
        sources=length_sources,
    )
    if not discovery_df.empty and not gene_lengths and args.filter_mode == 'stratified':
        print("Error: Could not load gene lengths from MAP directory or gene ensembl file", file=sys.stderr)
        sys.exit(1)

    discarded = read_discarded(args.cluster_file, discovery_df) if needs_clusters else set()
    print(f"Loaded {len(discovery_df)} CAAS positions across {discovery_df['Gene'].nunique()} genes", file=sys.stderr)

    units = build_units(discovery_df, discarded)
    removal = gene_removal(units, gene_lengths, args.filter_mode, strata=args.strata, k_far=args.k_far, k_near=args.k_near)
    removed_genes_df = pd.DataFrame(
        sorted((gene, cat) for (_lab, gene), cat in removal.items()),
        columns=['Gene', 'category'])

    filtered_df = apply_gene_filter(discovery_df, removed_genes_df)

    if args.remove_clusters:
        print("Applying cluster position removal (--remove-clusters enabled)...", file=sys.stderr)
        filtered_df = apply_cluster_filter(filtered_df, discarded)
    else:
        print("Cluster position removal: disabled (cluster positions retained unless removed by gene filter)", file=sys.stderr)

    print(f"Writing filtered discovery: {args.output}", file=sys.stderr)
    filtered_df.to_csv(args.output, sep='\t', index=False)
    if args.quality_universe and args.cluster_contract:
        import json
        from src.core.scoring_provenance import digest, save
        quality = json.loads(Path(args.quality_universe).read_text())
        cluster = json.loads(Path(args.cluster_contract).read_text())
        maps = cluster.pop('maps')
        cluster['remove'] = args.remove_clusters
        save(str(args.output)+'.scoring.json',dict(version=1,observed_sha256=digest(args.output),
             excluded_genes=quality.get('excluded_genes',[]),scientific=quality['scientific'],alignments=quality.get('alignments',{}),cluster=cluster,maps=maps,strict=quality.get('strict',False),
             gene_filter=dict(mode=args.filter_mode,strata=args.strata,k_far=args.k_far,k_near=args.k_near),lengths=gene_lengths,
             length_sources=length_sources))

    if args.summary:
        print(f"Writing removed genes summary: {args.summary}", file=sys.stderr)
        removed_genes_df.to_csv(args.summary, sep='\t', index=False)

    if args.gene_stats_output:
        stats = gene_unit_stats(units, gene_lengths, strata=args.strata, k_far=args.k_far, k_near=args.k_near)
        print(f"Writing gene statistics: {args.gene_stats_output}", file=sys.stderr)
        gene_stats_table(stats, removal, length_sources).to_csv(args.gene_stats_output, sep='\t', index=False)

    print("\n=== Post-Processing Filtering Summary ===", file=sys.stderr)
    print(f"Gene filter mode: {args.filter_mode}", file=sys.stderr)
    print(f"Cluster removal: {'enabled' if args.remove_clusters else 'disabled'}", file=sys.stderr)
    for category in ['Extreme', 'Dubious', 'Both']:
        count = int((removed_genes_df['category'] == category).sum())
        if count > 0:
            print(f"  {category}: {count}", file=sys.stderr)
    print(f"Total removed genes: {len(removed_genes_df)}", file=sys.stderr)
    pct_str = f"{100*len(filtered_df)/len(discovery_df):.2f}%" if len(discovery_df) > 0 else "N/A"
    print(f"Positions retained: {len(filtered_df)}/{len(discovery_df)} ({pct_str})", file=sys.stderr)


if __name__ == '__main__':
    main()
