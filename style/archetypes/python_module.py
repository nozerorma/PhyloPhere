#!/usr/bin/env python3
# filter_example.py — Filter discovery results by gene-level criteria.
# PhyloPhere | subworkflows/CT_POSTPROC/local/

"""
FilterExample: Applies dual filtering (density-percentile + IQR outliers) to
CAAS discovery output before downstream scoring and characterization.

Called by:  CT_POSTPROC Nextflow process (ctpp_clustfilter.nf → filter_caas_genes.py)
Inputs:     --discovery  TSV with [Trait, Gene, Position, ...] columns
            --lengths    TSV with [Gene, Length] for density normalization
            --output     Path for filtered output TSV
Outputs:    Filtered discovery TSV; log of excluded genes with reason
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import sys
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
import pandas as pd

# ── Local (only when part of a package — omit for true standalone scripts) ────
# from src.utils.io_utils import load_tsv


# ── Functions ─────────────────────────────────────────────────────────────────


def load_inputs(discovery_path: Path, lengths_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate input TSVs.

    Raises ValueError if required columns are missing rather than failing silently.
    This surfaces mismatches between pipeline versions early.
    """
    disc = pd.read_csv(discovery_path, sep="\t")
    lens = pd.read_csv(lengths_path, sep="\t")

    required_disc = {"Trait", "Gene", "Position"}
    missing = required_disc - set(disc.columns)
    if missing:
        raise ValueError(f"Discovery TSV missing columns: {missing}")

    return disc, lens


def filter_by_density(disc: pd.DataFrame, lens: pd.DataFrame, percentile: float = 0.99) -> pd.DataFrame:
    """Remove genes in the top `percentile` by CAAS/length density.

    High-density genes are likely alignment artifacts rather than true convergence
    signals; removing them reduces false positives in downstream scoring.
    """
    merged = disc.merge(lens, on="Gene", how="left")
    counts = disc.groupby("Gene").size().rename("n_caas")
    merged = merged.merge(counts, on="Gene", how="left")
    merged["density"] = merged["n_caas"] / merged["Length"].clip(lower=1)

    threshold = merged["density"].quantile(percentile)
    return merged[merged["density"] <= threshold].drop(columns=["density", "n_caas"])


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter CAAS discovery results by gene-level quality criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--discovery", required=True, type=Path,
                   help="Discovery TSV (output of ct_discovery)")
    p.add_argument("--lengths",   required=True, type=Path,
                   help="Gene length TSV [Gene, Length]")
    p.add_argument("--output",    required=True, type=Path,
                   help="Path for filtered output TSV")
    p.add_argument("--percentile", type=float, default=0.99,
                   help="Density percentile cutoff for extreme-gene removal (default: 0.99)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    disc, lens = load_inputs(args.discovery, args.lengths)
    filtered = filter_by_density(disc, lens, args.percentile)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"Retained {len(filtered)} / {len(disc)} positions after filtering.", file=sys.stderr)


if __name__ == "__main__":
    main()
