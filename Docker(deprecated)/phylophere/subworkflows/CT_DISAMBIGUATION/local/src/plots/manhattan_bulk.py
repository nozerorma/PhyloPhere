from pathlib import Path
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .plot_utils import chr_sort_key, compute_mean_evidence_score, safe_save

logger = logging.getLogger(__name__)


def plot_manhattan(df: pd.DataFrame, ensembl_df: pd.DataFrame, out: Path):
    """Generate Manhattan-style plot per chromosome using mean p-values."""
    if df.empty or ensembl_df.empty:
        logger.warning("No data available for manhattan plot")
        return
    if (
        "gene" not in df.columns
        or "gene" not in ensembl_df.columns
        or "chr" not in ensembl_df.columns
    ):
        logger.warning(
            "Required columns missing for manhattan plot (need gene and chr)"
        )
        return

    gene_stats = (
        df.groupby("gene")
        .apply(compute_mean_evidence_score)
        .reset_index(name="mean_pvalue")
        .dropna(subset=["mean_pvalue"])
    )

    # Merge with chromosome and optional position data
    merge_cols = ["gene", "chr"]
    if "start" in ensembl_df.columns and "end" in ensembl_df.columns:
        merge_cols.extend(["start", "end"])
        has_positions = True
    else:
        has_positions = False

    merged = gene_stats.merge(ensembl_df[merge_cols], on="gene", how="left")
    merged = merged.dropna(subset=["chr", "mean_pvalue"])
    if merged.empty:
        logger.warning(
            "No overlapping genes with chromosome data; skipping manhattan plot"
        )
        return

    merged["chr"] = merged["chr"].astype(str)

    # If we have position data, use within-chromosome positioning; otherwise use sequential
    if has_positions:
        merged["start"] = pd.to_numeric(merged["start"], errors="coerce")
        merged["end"] = pd.to_numeric(merged["end"], errors="coerce")
        # Position each gene at midpoint within chromosome
        merged["x_position"] = (merged["start"] + merged["end"]) / 2.0
        use_sequential = False
    else:
        use_sequential = True

    merged = merged.sort_values(by="chr", key=lambda col: col.map(chr_sort_key))

    # Use sequential indices (original behavior) if no position data, otherwise normalize per chromosome
    if use_sequential:
        merged["idx"] = range(len(merged))
    else:
        merged["idx"] = 0.0
        cumulative_offset = 0.0
        chr_width = 1.0  # give each chromosome the same visible width
        gap = 0.15  # small gap between chromosomes

        for chr_name in sorted(merged["chr"].unique(), key=chr_sort_key):
            chr_mask = merged["chr"] == chr_name
            chr_positions = merged.loc[chr_mask, "x_position"].values
            if chr_positions.size == 0:
                continue

            min_pos = chr_positions.min()
            max_pos = chr_positions.max()
            span = max(max_pos - min_pos, 1e-9)  # avoid divide-by-zero
            # Normalize start/end relative within chromosome so small chromosomes remain visible
            normalized_positions = (chr_positions - min_pos) / span
            merged.loc[chr_mask, "idx"] = (
                cumulative_offset + normalized_positions * chr_width
            )

            # Advance offset by fixed width plus gap so chromosomes don't overlap
            cumulative_offset += chr_width + gap

    merged["neglog_mean_p"] = merged["mean_pvalue"].apply(
        lambda v: -1 * np.log10(v) if v > 0 else np.nan
    )
    merged = merged.dropna(subset=["neglog_mean_p"])
    if merged.empty:
        logger.warning("Mean p-values are zero or invalid; skipping manhattan plot")
        return

    # Alternate colors by chromosome
    chromosomes = merged["chr"].unique()
    palette = ["#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974", "#64B5CD"]
    colors = {
        chr_val: palette[i % len(palette)] for i, chr_val in enumerate(chromosomes)
    }
    merged["color"] = merged["chr"].map(colors)

    fig, ax = plt.subplots(figsize=(16, 10))  # 3:2 aspect ratio
    ax.scatter(
        merged["idx"],
        merged["neglog_mean_p"],
        c=merged["color"],
        s=30,
        alpha=0.8,
        edgecolor="none",
    )
    ax.set_xlabel("Chromosome", fontweight="bold")
    ax.set_ylabel("-log10(mean_pvalue)", fontweight="bold")
    ax.set_title("Gene-level Manhattan plot (mean p-value)", fontweight="bold")

    # # Add Bonferroni threshold
    # num_genes = len(merged)
    # if num_genes > 0:
    #     bonferroni_pval = 0.05 / num_genes
    #     bonferroni_neglog = (
    #         -1 * np.log10(bonferroni_pval) if bonferroni_pval > 0 else np.inf
    #     )
    #     ax.axhline(
    #         y=bonferroni_neglog,
    #         color="red",
    #         linestyle="-",
    #         linewidth=2.5,
    #         alpha=0.8,
    #         label=f"Bonferroni threshold (α=0.05/{num_genes})",
    #         zorder=2,
    #     )
    #     ax.legend(loc="upper right", fontsize=10, framealpha=0.95)

    # Add chromosome delimiters (gray dashed vertical lines) and natural numbering (1,2,3...)
    chr_list = sorted(merged["chr"].unique(), key=chr_sort_key)
    chr_bounds = merged.groupby("chr")["idx"].agg(["min", "max"]).reindex(chr_list)
    chr_centers = chr_bounds.apply(lambda row: (row["min"] + row["max"]) / 2.0, axis=1)

    # Draw vertical delimiters between chromosome bands
    for i in range(len(chr_list) - 1):
        left_max = chr_bounds.iloc[i]["max"]
        right_min = chr_bounds.iloc[i + 1]["min"]
        boundary = (left_max + right_min) / 2.0
        ax.axvline(
            x=boundary, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, zorder=0
        )

    # Label with actual chromosome names at band centers
    x_tick_positions = [chr_centers[chr_val] for chr_val in chr_list]
    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels(chr_list, rotation=45, ha="right")

    # Add top 20 gene labels
    top_20 = merged.nlargest(20, "neglog_mean_p")
    for _, row in top_20.iterrows():
        ax.annotate(
            row["gene"],
            xy=(row["idx"], row["neglog_mean_p"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            alpha=0.8,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="yellow",
                alpha=0.3,
                edgecolor="none",
            ),
            arrowprops=dict(
                arrowstyle="->",
                connectionstyle="arc3,rad=0",
                color="gray",
                lw=0.8,
                alpha=0.5,
            ),
        )

    # Caption describing metric (moved down ~1cm = -0.22)
    # Composite evidence score for visuals
    # E = mean( -log10(E[p_value]), -log10(E[pvalue.boot]) )
    # i.e., mean of per-method mean -log10(p) across the group.
    if use_sequential:
        caption = (
            "evidence_score = ( mean( -log10(E[p_value]), -log10(E[pvalue.boot]) ) ); "
            "colors alternate by chromosome; sequential gene positioning"
        )
    else:
        caption = (
            "evidence_score = ( mean( -log10(E[p_value]), -log10(E[pvalue.boot]) ) ); "
            "colors alternate by chromosome; genes positioned relative within each chromosome"
        )

    ax.text(
        0.01, -0.22, caption, transform=ax.transAxes, fontsize=9, ha="left", va="top"
    )

    safe_save(fig, out)
