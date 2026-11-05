#!/usr/bin/env python3
# check_b0_removal.py — The genes removed from the null labeling b_0 must be the genes the observed filter removed.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Design invariant of the gene filter: b_0 is the real labeling replayed through the same tiles and merged with the same
rule and lengths as every null cycle, so the (gene, category) pairs it removes equal removed_genes_summary.tsv of the
observed filter. A difference means that the two sides evaluated different genes, lengths or counts, and the null is
then not comparable with the observed labeling.

Called by:  CAAS_REMOVAL_INVARIANT process (subworkflows/CT/caas_permulation.nf)
Inputs:     --b0        output directory of the CAAS_CORE_MERGE of b_0 (removed_units.tsv, units.parquet)
            --observed  removed_genes_summary.tsv of CAAS_FILTER_GENES (Gene, category)
Outputs:    --report    check_b0_removal.json; exit status 1 and a message naming the genes when the sets differ
"""

import argparse
import csv
import json
import sys
from pathlib import Path

OBSERVED_CYCLE = "b_0"


def read_observed(path):
    """removed_genes_summary.tsv to {gene: category}."""
    with open(path, newline="") as fh:
        return {row["Gene"]: row["category"] for row in csv.DictReader(fh, delimiter="\t")}


def read_b0(block_dir):
    """({gene: category}, categories known) of the genes removed from b_0; categories come from units.parquet."""
    block_dir = Path(block_dir)
    removed_path = block_dir / "removed_units.tsv"
    if not removed_path.is_file():
        return {}, True
    with open(removed_path, newline="") as fh:
        genes = {row["Gene"] for row in csv.DictReader(fh, delimiter="\t") if row["cycle"] == OBSERVED_CYCLE}
    categories = {g: None for g in genes}
    units = block_dir / "units.parquet"
    if genes and units.is_file():
        import pyarrow.parquet as pq
        for row in pq.read_table(units, columns=["gene", "cycle", "category"]).to_pylist():
            if row["cycle"] == OBSERVED_CYCLE and row["gene"] in categories:
                categories[row["gene"]] = row["category"]
    return categories, all(c is not None for c in categories.values())


def compare(b0, observed, with_category):
    """(only in the observed, only in b_0, genes removed by both with different categories)."""
    only_obs = sorted(set(observed) - set(b0))
    only_b0 = sorted(set(b0) - set(observed))
    differ = sorted(g for g in set(b0) & set(observed) if with_category and b0[g] != observed[g])
    return only_obs, only_b0, differ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--b0", required=True, type=Path)
    ap.add_argument("--observed", required=True, type=Path)
    ap.add_argument("--report", default="check_b0_removal.json", type=Path)
    args = ap.parse_args()
    observed = read_observed(args.observed)
    b0, with_category = read_b0(args.b0)
    only_obs, only_b0, differ = compare(b0, observed, with_category)
    args.report.write_text(json.dumps({"observed_removed": len(observed), "b0_removed": len(b0),
                                       "only_observed": only_obs, "only_b0": only_b0, "category_differs": differ,
                                       "categories_compared": with_category}, indent=1) + "\n")
    if only_obs or only_b0 or differ:
        def head(names):
            return ", ".join(names[:20]) + (f" ... (+{len(names) - 20})" if len(names) > 20 else "")
        print(f"ERROR: the gene removal of b_0 differs from the observed filter: observed {len(observed)} genes, b_0 {len(b0)} genes.\n"
              f"  removed only in the observed filter ({len(only_obs)}): {head(only_obs)}\n"
              f"  removed only in b_0 ({len(only_b0)}): {head(only_b0)}\n"
              f"  removed by both with another category ({len(differ)}): {head(differ)}\n"
              "The null would be filtered with a rule that does not reproduce the observed one.", file=sys.stderr)
        return 1
    print(f"b_0 removal equals the observed filter: {len(observed)} genes", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
