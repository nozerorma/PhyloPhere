#!/usr/bin/env python3
"""Validate an alignment and retain informative sequences present in a tree."""
import argparse
from Bio import SeqIO
from alignment_io import informative, read_alignment, tree_taxa


def filter_fasta(in_path, out_path, keep_taxa):
    records = read_alignment(in_path)
    kept = [r for r in records if informative(r) and r.id in keep_taxa]
    SeqIO.write(kept, out_path, "fasta")
    return len(kept), len(records) - len(kept)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    kept, _ = filter_fasta(args.fasta, args.output, tree_taxa(args.tree))
    if not kept:
        raise ValueError("No informative alignment sequences matched tree taxa")


if __name__ == "__main__":
    main()
