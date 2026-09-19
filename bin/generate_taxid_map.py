#!/usr/bin/env python3
"""
generate_taxid_map.py  —  Auto-generate a tax_id/species mapping from a
species tree's tip labels, using NCBI taxonomy for resolution.

Replaces the need for a user-supplied --tax_id file in the common case
(species names already match, or differ only by underscore/space).
Resolution is exact-match only: a tip label is looked up as its scientific
name (underscores -> spaces) against NCBI taxonomy. Anything that doesn't
resolve exactly is reported in --unresolved rather than guessed via fuzzy
or synonym matching, so the user can supply corrections deliberately.

Model: fam2ord.py (NEOPLASY_PRIMATES/Data/5.Phylogeny), which uses the same
ete3 NCBITaxa lookup pattern.

Output (--output): TSV with columns tax_id, species — the exact schema
required by subworkflows/TRAIT_ANALYSIS/local/src/phylo.R and
subworkflows/RERCONVERGE/local/rer_master_tree.R.
"""

import argparse
import csv
import re
import sys

import dendropy
from ete3 import NCBITaxa


def normalize_label(label: str) -> str:
    label = (label or "").strip().strip("'\"")
    return re.sub(r"\s+", "_", label)


def load_tip_labels(tree_path: str) -> list:
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    labels = [normalize_label(leaf.taxon.label) for leaf in tree.leaf_node_iter() if leaf.taxon]
    seen = set()
    unique = []
    for label in labels:
        if label not in seen:
            seen.add(label)
            unique.append(label)
    return unique


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick species tree")
    parser.add_argument("--output", required=True, help="Output TSV: tax_id, species")
    parser.add_argument("--unresolved", required=True,
                         help="Output TSV of tip labels that could not be resolved exactly")
    args = parser.parse_args()

    labels = load_tip_labels(args.tree)
    if not labels:
        print(f"Error: no tip labels found in {args.tree}", file=sys.stderr)
        sys.exit(1)

    ncbi = NCBITaxa()
    query_names = [label.replace("_", " ") for label in labels]
    name2taxid = ncbi.get_name_translator(query_names)

    resolved = []
    unresolved = []
    for label, query_name in zip(labels, query_names):
        hits = name2taxid.get(query_name)
        if hits:
            resolved.append((hits[0], label))
        else:
            unresolved.append(label)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["tax_id", "species"])
        for tax_id, species in resolved:
            writer.writerow([tax_id, species])

    with open(args.unresolved, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["species"])
        for species in unresolved:
            writer.writerow([species])

    print(f"Resolved {len(resolved)}/{len(labels)} tip labels against NCBI taxonomy.",
          file=sys.stderr)
    if unresolved:
        print(f"WARNING: {len(unresolved)} tip labels unresolved, see {args.unresolved}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
