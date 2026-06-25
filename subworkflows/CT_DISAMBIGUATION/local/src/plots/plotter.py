#!/usr/bin/env python3
"""
Bulk plotting utilities for CAAS aggregation outputs.

Generates visualizations from aggregated CAAS convergence results.
All plotting functions accept pre-loaded DataFrames and save results to disk.
"""

from pathlib import Path
import argparse
import logging
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt

from src.plots.gene_trees import plot_random_gene_trees

from src.plots.plot_utils import (
    load_df,
    find_tree_file,
)

logger = logging.getLogger(__name__)

import logging
import warnings

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
logging.getLogger("matplotlib.backends").setLevel(logging.WARNING)

warnings.filterwarnings("ignore", module="matplotlib")


def _generate_plot_suite(
    df: pd.DataFrame,
    output_dir: Path,
    n_gene_trees: int,
    asr_root: Path,
    node_dumps_root: Path,
    tip_details_root: Path,
):
    """Generate the full suite of plots."""
    sample_tree_plots_dir = output_dir / "sample_tree_plots"
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
    except Exception as exc:
        logger.debug("Could not emit per-gene debug info: %s", exc)

    plot_random_gene_trees(
        df,
        sample_tree_plots_dir,
        asr_root=asr_root,
        node_dumps_root=node_dumps_root,
        tip_details_root=tip_details_root,
        n=n_gene_trees,
    )


def generate_bulk_plots(
    caas_csv: Path,
    ensembl_csv: Optional[Path] = None,
    output_dir: Path = Path("plots"),
    n_gene_trees: int = 5,
    asr_root: Optional[Path] = None,
    node_dumps_root: Optional[Path] = None,
    tip_details_root: Optional[Path] = None,
):
    """Generate bulk plots from CAAS aggregation outputs (all results)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if asr_root is None:
        asr_root = output_dir.parent / "asr"
    if node_dumps_root is None:
        node_dumps_root = output_dir.parent / "diagnostics" / "node_dumps"
    if tip_details_root is None:
        tip_details_root = output_dir.parent / "diagnostics" / "tip_details"

    logger.info("Resolved ASR root: %s", asr_root)
    logger.info("Resolved node dumps root: %s", node_dumps_root)
    logger.info("Resolved tip details root: %s", tip_details_root)

    df = load_df(caas_csv)
    if df.empty:
        logger.warning("No CAAS data loaded")
        return

    logger.info("Plotting %s results", len(df))
    _generate_plot_suite(
        df,
        output_dir,
        n_gene_trees,
        asr_root,
        node_dumps_root,
        tip_details_root,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate bulk plots from CAAS aggregation"
    )
    parser.add_argument("caas_csv", type=Path, help="CAAS convergence CSV")
    parser.add_argument("--ensembl-csv", type=Path, help="Ensembl annotations CSV")
    parser.add_argument(
        "--output-dir", type=Path, default=Path("plots"), help="Output directory"
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
        args.n_gene_trees,
    )


if __name__ == "__main__":
    main()
