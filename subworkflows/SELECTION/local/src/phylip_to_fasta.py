#!/usr/bin/env python3
"""Normalize a validated alignment to FASTA, dropping uninformative sequences."""
import sys
from Bio import SeqIO
from alignment_io import informative, read_alignment


def phylip_to_fasta(input_path, output_path=None):
    records = [r for r in read_alignment(input_path) if informative(r)]
    if not records:
        raise ValueError(f"No informative sequences: {input_path}")
    SeqIO.write(records, output_path or sys.stdout, "fasta")


if __name__ == "__main__":
    if len(sys.argv) not in {2, 3}:
        sys.exit("Usage: phylip_to_fasta.py <input> [output.fa]")
    phylip_to_fasta(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else None)
