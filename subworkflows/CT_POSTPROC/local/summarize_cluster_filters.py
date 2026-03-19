#!/usr/bin/env python3
"""Summarize CT cluster-filter parameter sweeps."""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Build filter_summary.tsv and discarded_summary.tsv from cluster outputs."
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Directory containing *.filtered.*.tsv files",
    )
    parser.add_argument(
        "--summary-output",
        default="filter_summary.tsv",
        help="Output path for summary TSV",
    )
    parser.add_argument(
        "--discarded-output",
        default="discarded_summary.tsv",
        help="Output path for discarded-summary TSV",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    filter_files = sorted(Path(args.input_dir).glob("*.filtered.*.tsv"))
    if not filter_files:
        print("No filtered files found", file=sys.stderr)
        return 1

    df_sample = pd.read_csv(filter_files[0], sep="\t", nrows=5)
    has_caap = "CAAP_Group" in df_sample.columns

    print(f"Found {len(filter_files)} filtered files")
    print(f"CAAP mode: {has_caap}")

    summary_data = []
    discarded_data = []

    for file_path in filter_files:
        match = re.search(r"minlen(\d+)\.maxcaas(\d+)", file_path.stem)
        if not match:
            continue

        minlen = int(match.group(1))
        maxcaas = int(match.group(2))
        param_label = f"{minlen}/{maxcaas}"
        df = pd.read_csv(file_path, sep="\t")

        if has_caap and "CAAP_Group" in df.columns:
            discarded = df[df["ClusteringFlag"] == "Discarded"]
            group_counts = discarded.groupby("CAAP_Group").size().to_dict()
            row = {"Parameter": param_label, "Minlen": minlen, "Maxcaas": maxcaas}
            row.update({str(k): v for k, v in group_counts.items()})
            row["DiscardedCount"] = int(sum(group_counts.values()))
            summary_data.append(row)

            group_labels = [str(k) for k in group_counts.keys()]
            groups_str = ",".join(sorted(group_labels))
            discarded_data.append(
                {
                    "File": file_path.stem,
                    "DiscardedCount": row["DiscardedCount"],
                    "Minlen": minlen,
                    "Maxcaas": maxcaas,
                    "CAAP_Groups": groups_str if groups_str else "none",
                }
            )
        else:
            count = int((df["ClusteringFlag"] == "Discarded").sum())
            summary_data.append(
                {
                    "Parameter": param_label,
                    "Minlen": minlen,
                    "Maxcaas": maxcaas,
                    "DiscardedCount": count,
                }
            )
            discarded_data.append(
                {
                    "File": file_path.stem,
                    "DiscardedCount": count,
                    "Minlen": minlen,
                    "Maxcaas": maxcaas,
                    "CAAP_Groups": "N/A",
                }
            )

    if has_caap:
        summary_df = pd.DataFrame(summary_data).fillna(0)
        summary_df = summary_df.sort_values(["Minlen", "Maxcaas"])
        protected = ["Parameter", "Minlen", "Maxcaas", "DiscardedCount"]
        group_cols = sorted(col for col in summary_df.columns if col not in protected)
        summary_df = summary_df[protected + group_cols]
        for col in ["DiscardedCount", *group_cols]:
            summary_df[col] = summary_df[col].astype(int)
    else:
        summary_df = pd.DataFrame(summary_data).sort_values(["Minlen", "Maxcaas"])
        summary_df["DiscardedCount"] = summary_df["DiscardedCount"].astype(int)

    discarded_df = pd.DataFrame(discarded_data).sort_values(["Minlen", "Maxcaas"])
    discarded_df["DiscardedCount"] = discarded_df["DiscardedCount"].astype(int)

    summary_df.to_csv(args.summary_output, sep="\t", index=False)
    discarded_df.to_csv(args.discarded_output, sep="\t", index=False)

    print(f"✓ {args.summary_output}: {len(summary_df)} parameter combinations")
    print(f"✓ {args.discarded_output}: {len(discarded_df)} entries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
