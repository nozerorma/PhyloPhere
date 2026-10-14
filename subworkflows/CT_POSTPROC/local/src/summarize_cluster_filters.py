#!/usr/bin/env python3
# summarize_cluster_filters.py — Count the positions discarded by each (minlen, maxcaas) cluster filter.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Summary of the cluster-filter parameter sweep: the number of Discarded positions of every
*.filtered.minlen<L>.maxcaas<C>.tsv file.

Called by:  CT_FILTER_SUMMARY process (ctpp_clustfilter.nf)
Inputs:     --input-dir  directory with the *.filtered.*.tsv files of filter_caas_clusters-param.py
Outputs:    --summary-output  one row per parameter pair (Parameter, Minlen, Maxcaas, DiscardedCount)
            --discarded-output  one row per parameter pair (File, DiscardedCount, Minlen, Maxcaas)
"""

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

    print(f"Found {len(filter_files)} filtered files")

    summary_data = []

    for file_path in filter_files:
        match = re.search(r"minlen(\d+)\.maxcaas(\d+)", file_path.stem)
        if not match:
            continue

        minlen = int(match.group(1))
        maxcaas = int(match.group(2))
        df = pd.read_csv(file_path, sep="\t")
        summary_data.append(
            {
                "Parameter": f"{minlen}/{maxcaas}",
                "Minlen": minlen,
                "Maxcaas": maxcaas,
                "DiscardedCount": int((df["clustering_flag"] == "Discarded").sum()),
            }
        )

    summary_df = pd.DataFrame(summary_data).sort_values(["Minlen", "Maxcaas"])
    discarded_df = summary_df[["Parameter", "DiscardedCount", "Minlen", "Maxcaas"]].rename(columns={"Parameter": "File"})

    summary_df.to_csv(args.summary_output, sep="\t", index=False)
    discarded_df.to_csv(args.discarded_output, sep="\t", index=False)

    print(f"✓ {args.summary_output}: {len(summary_df)} parameter combinations")
    print(f"✓ {args.discarded_output}: {len(discarded_df)} entries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
