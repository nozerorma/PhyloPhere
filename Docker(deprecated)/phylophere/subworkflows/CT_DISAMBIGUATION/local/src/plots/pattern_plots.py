from pathlib import Path
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .plot_utils import get_pattern_abbreviation, get_pattern_colors, safe_save

logger = logging.getLogger(__name__)


def _plot_pattern_stack(
    df: pd.DataFrame, out: Path, title: str, theme_color: str = "#16a085"
):
    counts = df["pattern_type"].value_counts().sort_values(ascending=False)
    pattern_order = list(counts.index)
    if not pattern_order:
        logger.warning("No pattern_type values available for %s", out)
        return

    x_pos = np.arange(len(pattern_order))
    labels = [get_pattern_abbreviation(p) for p in pattern_order]

    raw_counts = [counts.get(p, 0) for p in pattern_order]
    sig_counts = []
    asr_cons_counts = []

    for pattern in pattern_order:
        subset = df[df["pattern_type"] == pattern]
        sig = (
            subset.get("is_significant", pd.Series(True, index=subset.index))
            .fillna(False)
            .astype(bool)
        )
        sig_counts.append(sig.sum())

        asr_c = (
            subset.get("asr_is_conserved", pd.Series(False, index=subset.index))
            .fillna(False)
            .astype(bool)
            .sum()
        )
        asr_cons_counts.append(asr_c)

    fig, axes = plt.subplots(3, 1, figsize=(16, 10))
    fig.suptitle(title, fontsize=13, fontweight="bold")

    axes[0].bar(
        x_pos,
        raw_counts,
        color=[get_pattern_colors(p) for p in pattern_order],
        edgecolor="black",
        linewidth=1.0,
    )
    axes[0].set_title("Row 1: Raw pattern counts", loc="left", fontweight="bold")
    axes[0].set_ylabel("Count")
    axes[0].set_xticks([])
    axes[0].grid(axis="y", alpha=0.3)

    non_sig = [max(raw_counts[i] - sig_counts[i], 0) for i in range(len(pattern_order))]
    axes[1].bar(
        x_pos,
        non_sig,
        color="#ecf0f1",
        edgecolor="black",
        linewidth=1.0,
        label="Non-significant",
    )
    axes[1].bar(
        x_pos,
        sig_counts,
        bottom=non_sig,
        color=theme_color,
        edgecolor="black",
        linewidth=1.0,
        label="Significant",
    )
    axes[1].set_title("Row 2: Significance overlay", loc="left", fontweight="bold")
    axes[1].set_ylabel("Count")
    axes[1].legend(loc="upper right", fontsize=9)
    axes[1].set_xticks([])
    axes[1].grid(axis="y", alpha=0.3)

    other = [max(raw_counts[i] - sig_counts[i], 0) for i in range(len(pattern_order))]
    axes[2].bar(
        x_pos, other, color="#3498db", edgecolor="black", linewidth=1.0, label="Other"
    )
    axes[2].bar(
        x_pos,
        sig_counts,
        bottom=other,
        color="#f39c12",
        edgecolor="black",
        linewidth=1.0,
        label="Significant",
    )
    axes[2].bar(
        x_pos,
        asr_cons_counts,
        bottom=[other[i] + sig_counts[i] for i in range(len(pattern_order))],
        color="#27ae60",
        edgecolor="black",
        linewidth=1.0,
        label="ASR conserved",
    )
    axes[2].set_title(
        "Row 3: Significant + ASR-conserved overlay", loc="left", fontweight="bold"
    )
    axes[2].set_ylabel("Count")
    axes[2].set_xlabel("Pattern type")
    axes[2].legend(loc="upper right", fontsize=9)
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(labels, rotation=0, ha="center", fontweight="bold")
    axes[2].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    safe_save(fig, out)


def plot_convergence_patterns(df: pd.DataFrame, out: Path):
    _plot_pattern_stack(df, out, "Pattern classification overview")


def plot_faceted_patterns(df: pd.DataFrame, out: Path):
    facets = {
        "top": (df["change_side"].isin(["top", "both"]), "TOP/BOTH changes", "#cc0000"),
        "bottom": (
            df["change_side"].isin(["bottom", "both"]),
            "BOTTOM/BOTH changes",
            "#00aa00",
        ),
        "both": (df["change_side"] == "both", "BOTH-side changes", "#7d3c98"),
    }
    for key, (mask, title, color) in facets.items():
        subset = df[mask].copy()
        if subset.empty:
            continue
        _plot_pattern_stack(
            subset,
            out.parent / f"{out.stem}_facet_{key}.png",
            f"Faceted patterns: {title}",
            theme_color=color,
        )


def plot_trait_high_low_patterns(df: pd.DataFrame, out: Path, trait_type: str):
    _plot_pattern_stack(
        df,
        out,
        f"Trait-specific patterns ({trait_type})",
        theme_color="#cc0000" if trait_type == "high" else "#00aa00",
    )
