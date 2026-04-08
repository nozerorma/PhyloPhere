#!/usr/bin/env python3
"""
Filter a FASTA alignment down to taxa present in a reference tree.
"""

import argparse
import sys

import dendropy


def tree_taxa(path: str) -> set:
    with open(path) as fh:
        newick = fh.read().strip()
    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    return {leaf.taxon.label for leaf in tree.leaf_node_iter()}


def filter_fasta(in_path: str, out_path: str, keep_taxa: set) -> tuple[int, int]:
    kept = 0
    dropped = 0
    include = False
    with open(in_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith(">"):
                taxon = line[1:].strip().split()[0]
                include = taxon in keep_taxa
                if include:
                    kept += 1
                else:
                    dropped += 1
            if include:
                fh_out.write(line)
    return kept, dropped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter a FASTA alignment to the taxa present in a tree"
    )
    parser.add_argument("--tree", required=True, help="Reference Newick tree")
    parser.add_argument("--fasta", required=True, help="Input FASTA alignment")
    parser.add_argument("--output", required=True, help="Filtered FASTA output path")
    args = parser.parse_args()

    keep_taxa = tree_taxa(args.tree)
    kept, dropped = filter_fasta(args.fasta, args.output, keep_taxa)

    if kept == 0:
        sys.exit(
            "ERROR filter_fasta_to_tree: no alignment sequences matched taxa in the reference tree."
        )

    sys.stderr.write(
        f"INFO filter_fasta_to_tree: kept {kept} sequence(s), dropped {dropped} "
        f"sequence(s) absent from the reference tree.\n"
    )


if __name__ == "__main__":
    main()
