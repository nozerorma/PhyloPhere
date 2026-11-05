#!/usr/bin/env python3
# check_b0_removal.py — The genes removed from the null labeling b_0 must be the genes the observed filter removed.
# PhyloPhere | subworkflows/CT_POSTPROC/local/src/

"""
Design invariant of the gene filter: b_0 is the real labeling replayed through the same tiles and merged with the same
rule and lengths as every null cycle, so the (gene, category) pairs it removes equal removed_genes_summary.tsv of the
observed filter. A difference means that the two sides evaluated different genes, lengths or counts, and the null is
then not comparable with the observed labeling.

The settings come first: the null records the cluster parameters, the gene-filter settings and the strict flag it was
built with (counts.json of b_0), and SCORING requires the observed contract to equal them. With --observed-contract the
check names the settings that differ (for example, gene-removal switches changed after the null was merged) before it
compares the genes.

The b_0 directory is the output of this run's CAAS_CORE_MERGE, or the b0/ directory of a reused null
(<caas_permulation>/b0, <outdir>/strict/caas_permulation/b0 for the strict null), so a pass that merges nothing is checked too.

Called by:  CAAS_REMOVAL_INVARIANT process (subworkflows/CT/caas_permulation.nf)
Inputs:     --b0        directory of the merged b_0 (removed_units.tsv, units.parquet, counts.json)
            --observed  removed_genes_summary.tsv of CAAS_FILTER_GENES (Gene, category)
            --observed-contract  filtered_discovery.tsv.scoring.json of CAAS_FILTER_GENES (optional)
            --settings-only      compare the settings and not the genes: for a run whose observed table is a given
                                 discovery.tab, which is not the labeling the null replayed
Outputs:    --report    check_b0_removal.json; exit status 1 and a message naming the settings or genes that differ
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


SETTINGS = ("cluster", "gene_filter", "strict")


def settings_differences(b0_dir, contract_path):
    """Messages for the settings recorded in b_0 (counts.json) that differ from the observed contract; empty when equal.

    A b_0 without scoring provenance (tiles written without it) has nothing to compare and gives no message.
    """
    counts_path = Path(b0_dir) / "counts.json"
    if not counts_path.is_file():
        return []
    provenance = json.loads(counts_path.read_text()).get("scoring_provenance")
    if not isinstance(provenance, dict):
        return []
    contract = json.loads(Path(contract_path).read_text())
    messages = []
    for key in SETTINGS:
        null_value, observed_value = provenance.get(key), contract.get(key)
        if null_value != observed_value:
            messages.append(f"{key}: the null was built with {json.dumps(null_value, sort_keys=True)} but the observed filter used "
                            f"{json.dumps(observed_value, sort_keys=True)}")
    return messages


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
    ap.add_argument("--observed-contract", type=Path, default=None)
    ap.add_argument("--settings-only", action="store_true")
    ap.add_argument("--report", default="check_b0_removal.json", type=Path)
    args = ap.parse_args()
    differences = settings_differences(args.b0, args.observed_contract) if args.observed_contract else []
    if differences:
        args.report.write_text(json.dumps({"settings_differ": differences}, indent=1) + "\n")
        print("ERROR: the null and the observed filter were built with different settings:\n  " + "\n  ".join(differences) +
              "\nChanging a gene-removal setting needs the merge of the null to be rerun with it (-resume recomputes only the merges); "
              "changing a cluster parameter recomputes the tiles.", file=sys.stderr)
        return 1
    if args.settings_only:
        args.report.write_text(json.dumps({"settings_differ": [], "genes_compared": False}, indent=1) + "\n")
        print("the settings of the null equal those of the observed filter (genes not compared)", file=sys.stderr)
        return 0
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
