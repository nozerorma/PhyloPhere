#!/usr/bin/env python3
"""Normalize ct_disambiguation output for CT post-processing."""

import argparse
import re
import sys

import pandas as pd


def _normalize_schema(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "gene" in df.columns:
        rename_map["gene"] = "Gene"
    if "msa_pos" in df.columns:
        rename_map["msa_pos"] = "Position"
    if "caap_group" in df.columns:
        rename_map["caap_group"] = "CAAP_Group"
    if "pattern_type" in df.columns:
        rename_map["pattern_type"] = "Pattern"
    if "pvalue" in df.columns:
        rename_map["pvalue"] = "Pvalue"
    if rename_map:
        df = df.rename(columns=rename_map)

    if "Trait" not in df.columns:
        df["Trait"] = "post_disambiguation"
    if "CAAP_Group" not in df.columns:
        df["CAAP_Group"] = "US"
    if "is_conserved_meta" not in df.columns:
        df["is_conserved_meta"] = False
    if "asr_is_conserved" not in df.columns:
        df["asr_is_conserved"] = False

    return df


def _collect_removed_rows(cleaned: pd.DataFrame, mrca_threshold: float):
    removed_frames = []

    mrca_cols = [
        col for col in cleaned.columns if re.fullmatch(r"mrca_\d+_posterior", str(col))
    ]
    removed_low_mrca = cleaned.iloc[0:0].copy()
    if mrca_cols:
        for col in mrca_cols:
            cleaned[col] = pd.to_numeric(cleaned[col], errors="coerce")
        mask_low_mrca = cleaned[mrca_cols].lt(mrca_threshold).any(axis=1, skipna=True)
        removed_low_mrca = cleaned.loc[mask_low_mrca].copy()
        if not removed_low_mrca.empty:
            removed_low_mrca["RemovalReason"] = (
                f"mrca_posterior_below_{mrca_threshold:g}"
            )
            removed_frames.append(removed_low_mrca)
        cleaned = cleaned.loc[~mask_low_mrca].copy()

    return cleaned, removed_frames, removed_low_mrca, mrca_cols


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare ct_disambiguation master CSV for CT post-processing."
    )
    parser.add_argument("--input", required=True, help="Disambiguation master CSV/TSV")
    parser.add_argument(
        "--mrca-threshold",
        required=True,
        type=float,
        help="Canonical posterior threshold for mrca_*_posterior filtering",
    )
    parser.add_argument(
        "--output",
        default="postproc_disambiguation_input.tsv",
        help="Normalized output TSV",
    )
    parser.add_argument(
        "--removed-output",
        default="removed_patterns_precluster.tsv",
        help="Precluster removal output TSV",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep=None, engine="python")
    df = _normalize_schema(df)

    for col in ("Gene", "Position"):
        if col not in df.columns:
            raise ValueError(f"Required column missing after normalization: {col}")

    cleaned, removed_frames, removed_low_mrca, mrca_cols = (
        _collect_removed_rows(df.copy(), args.mrca_threshold)
    )

    cleaned["Position"] = pd.to_numeric(cleaned["Position"], errors="raise").astype(int)

    if removed_frames:
        removed = pd.concat(removed_frames, ignore_index=True)
        if "Position" in removed.columns:
            removed["Position"] = pd.to_numeric(removed["Position"], errors="coerce")
    else:
        removed = pd.DataFrame(columns=[*df.columns, "RemovalReason"])

    cleaned.to_csv(args.output, sep="\t", index=False)
    removed.to_csv(args.removed_output, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    if mrca_cols:
        print(
            f"Removed low MRCA posterior (< {args.mrca_threshold}): {len(removed_low_mrca)}"
        )
        print(f"MRCA posterior columns used: {', '.join(mrca_cols)}")
    else:
        print("No mrca_*_posterior columns found; MRCA posterior pruning skipped")
    print(f"Rows kept for postproc filtering: {len(cleaned)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
