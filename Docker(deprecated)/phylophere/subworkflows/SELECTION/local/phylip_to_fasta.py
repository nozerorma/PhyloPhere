#!/usr/bin/env python3
"""
phylip_to_fasta.py
──────────────────
Convert a relaxed-PHYLIP protein alignment to FASTA format.

Usage
-----
    python phylip_to_fasta.py <input.phy> [output.fa]

If output path is omitted the FASTA is written to stdout.

Notes
-----
* Uses BioPython's AlignIO which handles both strict and relaxed PHYLIP.
* Accepts both interleaved and sequential PHYLIP files.
* Gap-only sequences are silently dropped to avoid FADE/MoleRate failures.
"""

import sys
import os
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def phylip_to_fasta(input_path: str, output_path: str | None = None) -> None:
    """Read a PHYLIP alignment and write it as FASTA."""

    # Try relaxed first, then sequential as fallback
    for fmt in ("phylip-relaxed", "phylip-sequential", "phylip"):
        try:
            alignment = AlignIO.read(input_path, fmt)
            break
        except Exception:
            alignment = None

    if alignment is None:
        sys.exit(f"ERROR phylip_to_fasta: could not parse '{input_path}' as any PHYLIP variant")

    records = []
    for rec in alignment:
        # Strip internal whitespace from sequence (some PHYLIP writers leave spaces)
        clean_seq = str(rec.seq).replace(" ", "").replace("\t", "")
        # Drop gap-only sequences
        non_gap = clean_seq.replace("-", "").replace("?", "").replace("X", "").replace("x", "")
        if not non_gap:
            sys.stderr.write(f"WARNING: dropping gap-only sequence '{rec.id}'\n")
            continue
        records.append(SeqRecord(Seq(clean_seq), id=rec.id, description=""))

    if not records:
        sys.exit(f"ERROR phylip_to_fasta: no valid sequences in '{input_path}'")

    if output_path:
        SeqIO.write(records, output_path, "fasta")
    else:
        for rec in records:
            print(f">{rec.id}")
            print(str(rec.seq))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: phylip_to_fasta.py <input.phy> [output.fa]")
    inp  = sys.argv[1]
    outp = sys.argv[2] if len(sys.argv) > 2 else None
    phylip_to_fasta(inp, outp)
