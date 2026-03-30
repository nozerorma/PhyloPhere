#!/usr/bin/env python3
"""
Bulk plotting utilities for CAAS aggregation outputs.

This module generates comprehensive visualization sets from aggregated CAAS convergence
results. All plotting functions accept pre-loaded DataFrames
(pure, testable functions) and save results to disk.

**Design principles:**
- Pure functions accepting DataFrames as input (no hidden I/O)
- Deterministic ordering and reproducible plots
- Comprehensive logging for debugging and traceability

**Generated plots:**

1. **Pattern-based plots:**
   - convergence_patterns.png: 3-row structure (raw, significance, significant + ASR-conserved)
   - convergence_patterns_faceted_*.png: Faceted by change_side (top/bottom/both)
   - convergence_patterns_trait_*.png: Filtered by trait_value (high/low)

2. **Supporting plots:**
   - trait_*_analysis.png: Original 2-panel trait analysis

Author
------
Miguel Ramon Alonso
Evolutionary Genomics Lab - IBE-UPF

Date: 2025-12-07
Updated: 2026-02-12 (Simplified plots for convergence-only pipeline)
"""

from pathlib import Path
import argparse
import logging
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt

from src.plots.gene_trees_bulk import plot_random_gene_trees
from src.plots.manhattan_bulk import plot_manhattan
from src.plots.pattern_plots import (
    plot_convergence_patterns,
    plot_faceted_patterns,
    plot_trait_high_low_patterns,
)
from src.plots.plot_utils import (
    load_df,
    safe_save,
    get_pattern_abbreviation,
    find_tree_file,
)

logger = logging.getLogger(__name__)

import logging
import warnings

# Suppress matplotlib log spam
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("matplotlib.backends").setLevel(logging.WARNING)

# Optional: suppress matplotlib-related warnings
warnings.filterwarnings("ignore", module="matplotlib")


def _generate_plot_suite(
    df: pd.DataFrame,
    ensembl_df: pd.DataFrame,
    output_dir: Path,
    n_gene_trees: int,
    asr_root: Path,
    node_dumps_root: Path,
    tip_details_root: Path,
):
    """Generate the full suite of bulk plots."""
    # Manhattan plot
    manhattan_out = output_dir / "manhattan_mean_pvalue.png"
    plot_manhattan(df, ensembl_df, manhattan_out)

    # Convergence patterns
    patterns_out = output_dir / "convergence_patterns.png"
    plot_convergence_patterns(df, patterns_out)

    # Faceted patterns
    plot_faceted_patterns(df, output_dir / "convergence_patterns_faceted")

    # Trait patterns
    if "trait_value" in df.columns:
        for trait in ["high", "low"]:
            trait_df = df[df["trait_value"] == trait]
            if not trait_df.empty:
                trait_out = output_dir / f"convergence_patterns_trait_{trait}.png"
                plot_trait_high_low_patterns(trait_df, trait_out, trait)

    # Trait analysis
    if "change_side" in df.columns:
        for focus in ["top", "bottom"]:
            trait_analysis_out = output_dir / f"trait_{focus}_analysis.png"
            plot_trait_analysis(df, trait_analysis_out, focus)

    # Gene trees
    gene_trees_dir = output_dir / "gene_tree_samples"
    # Provide additional per-gene debug logging so users can see why gene trees
    try:
        unique_genes = []
        if "gene" in df.columns:
            unique_genes = list(pd.unique(df["gene"]))
        sample_genes = unique_genes[: max(0, min(len(unique_genes), n_gene_trees * 3))]
        if sample_genes:
            logger.debug("Sample genes for tree plotting: %s", sample_genes)
            for g in sample_genes:
                try:
                    tree_path = find_tree_file(str(g), asr_root)
                except Exception as e:
                    tree_path = None
                    logger.debug("find_tree_file raised for %s: %s", g, e)

                post_path = node_dumps_root / f"{str(g).lower()}_posteriors.jsonl"
                logger.debug(
                    "Gene %s: tree=%s posteriors=%s",
                    g,
                    str(tree_path) if tree_path else "(none)",
                    "exists" if post_path.exists() else "(none)",
                )

                tip_path = tip_details_root / f"{str(g).lower()}_tip_details.jsonl"
                logger.debug(
                    "Gene %s: tip details=%s",
                    g,
                    "exists" if tip_path.exists() else "(none)",
                )
    except Exception as exc:  # defensive - don't stop plotting on logging failure
        logger.debug("Could not emit per-gene debug info: %s", exc)

    plot_random_gene_trees(
        df,
        gene_trees_dir,
        asr_root=asr_root,
        node_dumps_root=node_dumps_root,
        tip_details_root=tip_details_root,
        n=n_gene_trees,
    )


def generate_bulk_plots(
    caas_csv: Path,
    ensembl_csv: Optional[Path] = None,
    output_dir: Path = Path("plots_bulk"),
    include_non_significant: bool = False,
    n_gene_trees: int = 5,
    asr_root: Optional[Path] = None,
    node_dumps_root: Optional[Path] = None,
    tip_details_root: Optional[Path] = None,
):
    """Generate comprehensive bulk plots from CAAS aggregation outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve defaults for ASR and node dump roots (keep legacy locations)
    if asr_root is None:
        asr_root = output_dir.parent / "asr"
    if node_dumps_root is None:
        node_dumps_root = output_dir.parent / "diagnostics" / "node_dumps"
    if tip_details_root is None:
        tip_details_root = output_dir.parent / "diagnostics" / "tip_details"

    logger.info("Resolved ASR root: %s", asr_root)
    logger.info("Resolved node dumps root: %s", node_dumps_root)
    logger.info("Resolved tip details root: %s", tip_details_root)

    # Load data
    df_full = load_df(caas_csv)
    if df_full.empty:
        logger.warning("No CAAS data loaded")
        return

    ensembl_df = load_df(ensembl_csv) if ensembl_csv else pd.DataFrame()

    # Prepare significant subset (do not overwrite df_full so we can optionally
    # still generate 'all' outputs when requested).
    if "is_significant" in df_full.columns:
        df_significant = df_full[
            df_full["is_significant"].fillna(False).astype(bool)
        ].copy()
        logger.info(
            "Significant filter: %s/%s results", len(df_significant), len(df_full)
        )
    else:
        df_significant = df_full.copy()
        logger.warning(
            "No 'is_significant' column found; treating all results as significant"
        )

    # Generate plots for significant results in 'significant' subdirectory
    significant_dir = output_dir / "significant"
    significant_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Generating plots for significant results in %s...", significant_dir)
    _generate_plot_suite(
        df_significant,  # type : ignore
        ensembl_df,
        significant_dir,
        n_gene_trees,
        asr_root,
        node_dumps_root,
        tip_details_root,
    )

    # If requested, also generate plots for ALL results in 'all' subdirectory
    if include_non_significant:
        all_dir = output_dir / "all"
        all_dir.mkdir(parents=True, exist_ok=True)
        logger.info(
            "Generating plots for all results (including non-significant) in %s...",
            all_dir,
        )
        _generate_plot_suite(
            df_full,
            ensembl_df,
            all_dir,
            n_gene_trees,
            asr_root,
            node_dumps_root,
            tip_details_root,
        )


def plot_trait_analysis(df: pd.DataFrame, out: Path, focus: str):
    """Plot trait-specific pattern distributions for a focus side."""

    subset = df[df["change_side"].isin([focus, "both"])]
    if subset.empty:
        logger.warning("No data for trait_%s_analysis", focus)
        return

    subset = subset.copy()

    counts = subset["pattern_type"].value_counts().sort_values(ascending=False)

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 10))  # 3:2 aspect ratio

    # Left: raw pattern distribution
    patterns = list(counts.index)
    x = range(len(patterns))

    axes[0].bar(x, counts.values, color="teal", label="All")

    axes[0].set_xticks(list(x))
    pattern_labels = [get_pattern_abbreviation(p) for p in patterns]
    axes[0].set_xticklabels(pattern_labels, rotation=0, ha="center", fontweight="bold")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"Trait-{focus.upper()} Pattern Distribution (All)")
    axes[0].set_xlabel("Pattern Type")

    # Right: significant subset by pattern
    sig_mask = (
        subset.get("is_significant", pd.Series(False, index=subset.index))
        .fillna(False)
        .astype(bool)
    )
    sig_counts = (
        subset[sig_mask].groupby("pattern_type").size().reindex(patterns, fill_value=0)
    )
    axes[1].bar(x, sig_counts.values, color="#16a085")
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(pattern_labels, rotation=0, ha="center", fontweight="bold")
    axes[1].set_title(f"Trait-{focus.upper()} Significant Subset")
    axes[1].set_ylabel("Count")
    axes[1].set_xlabel("Pattern Type")

    safe_save(fig, out)


def main():
    parser = argparse.ArgumentParser(
        description="Generate bulk plots from CAAS aggregation"
    )
    parser.add_argument("caas_csv", type=Path, help="CAAS convergence CSV")
    parser.add_argument("--ensembl-csv", type=Path, help="Ensembl annotations CSV")
    parser.add_argument(
        "--output-dir", type=Path, default=Path("plots_bulk"), help="Output directory"
    )
    parser.add_argument(
        "--include-non-significant",
        action="store_true",
        help="Include non-significant results",
    )
    parser.add_argument(
        "--n-gene-trees", type=int, default=5, help="Number of gene trees to plot"
    )

    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )

    generate_bulk_plots(
        args.caas_csv,
        args.ensembl_csv,
        args.output_dir,
        args.include_non_significant,
        args.n_gene_trees,
    )


if __name__ == "__main__":
    main()
