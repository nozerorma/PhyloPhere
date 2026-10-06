#!/usr/bin/env python3
# prepare_postproc_input.py — Normalize the ct_disambiguation master table for CT post-processing.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Input preparation of CT post-processing: standardizes the key columns of the disambiguation
master table and, when the alignments and contrast species are given, adds the extant-species
residue tally (residue_descriptors.py).

Called by:  CAAS_PREPARE_POSTPROC_INPUT process (ctpp_clustfilter.nf)
Inputs:     --input  disambiguation master CSV/TSV (the separator is detected from the header)
            --alignment, --alignment-format, --fg-species, --bg-species  optional, for the tally
Outputs:    --output  normalized TSV (Gene, Position, caap_group, ...)
            --removed-output  TSV of precluster removals (header only: no row is removed)
"""

import argparse
import os
import sys

import pandas as pd

from residue_descriptors import add_species_tally


def _normalize_schema(df: pd.DataFrame) -> pd.DataFrame:
    # Only the structural keys are renamed (gene/msa_pos → Gene/Position, the names the
    # downstream steps read). The other columns (caap_group, convergence_type, pvalue, caas,
    # amino_encoded) keep their lowercase names.
    rename_map = {}
    if "gene" in df.columns:
        rename_map["gene"] = "Gene"
    if "msa_pos" in df.columns:
        rename_map["msa_pos"] = "Position"
    if rename_map:
        df = df.rename(columns=rename_map)

    if "caap_group" not in df.columns:
        df["caap_group"] = "US"

    return df


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare ct_disambiguation master CSV for CT post-processing."
    )
    parser.add_argument("--input", required=True, help="Disambiguation master CSV/TSV")
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
    # Optional extant-species residue tally (top/bottom_species_residues, n_top/bottom_species).
    # It needs the alignment directory and both contrast species lists; without them the
    # columns stay empty.
    parser.add_argument("--alignment", default=None,
                        help="Alignment directory (flat, gene = basename up to first '.')")
    parser.add_argument("--alignment-format", default="fasta",
                        help="Bio.AlignIO format for --alignment (default: fasta)")
    parser.add_argument("--fg-species", default=None,
                        help="top_species.txt (foreground contrast species, one per line)")
    parser.add_argument("--bg-species", default=None,
                        help="bottom_species.txt (background contrast species, one per line)")
    args = parser.parse_args()

    # keep_default_na=False + na_values=[""]: the disambiguation master has categorical
    # amino-acid columns (caas, amino_encoded, mrca_*_aa) that can equal NA-sentinel
    # strings; only an empty cell is missing here.
    # The C engine with float_precision="round_trip" keeps every bit of the float columns
    # (asr_path_score feeds a gene score that compares values exactly); the python engine
    # that autodetects the separator does not. The separator is read from the header line.
    with open(args.input, newline="") as fh:
        header = fh.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    df = pd.read_csv(args.input, sep=sep, engine="c", float_precision="round_trip",
                     keep_default_na=False, na_values=["", "nan", "NaN"])
    df = _normalize_schema(df)

    for col in ("Gene", "Position"):
        if col not in df.columns:
            raise ValueError(f"Required column missing after normalization: {col}")

    cleaned = df.copy()
    cleaned["Position"] = pd.to_numeric(cleaned["Position"], errors="raise").astype(int)

    # Extant-species residue tally: a no-op when the alignment or the species lists are not
    # supplied (e.g. a standalone --disambiguation_input run).
    def _opt(p):
        return p if p and os.path.exists(p) else None
    cleaned = add_species_tally(
        cleaned,
        _opt(args.alignment),
        _opt(args.fg_species),
        _opt(args.bg_species),
        ali_format=args.alignment_format,
    )

    removed = pd.DataFrame(columns=[*df.columns, "removal_reason"])

    cleaned.to_csv(args.output, sep="\t", index=False)
    removed.to_csv(args.removed_output, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print("Precluster MRCA posterior pruning retired; all input rows kept for postproc filtering")
    print(f"Rows kept for postproc filtering: {len(cleaned)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
