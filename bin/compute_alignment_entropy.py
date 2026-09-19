#!/usr/bin/env python3
"""
compute_alignment_entropy.py  —  Batch driver auto-generating the per-gene
Valdar variability files CT_ACCUMULATION's --accumulation_entropy_dir expects,
directly from the alignment directory (no cds2prot/ortholog_characterizator
run required).

bin/compute_variability.py in this same directory is a byte-for-byte copy of
ortholog_characterizator/subworkflows/variability/local/compute_variability.py
(the user's own reference implementation) — fully verbatim, not
reimplemented, so the numbers this produces match what that pipeline would
have produced for the same alignment. This driver only adds the batching
(the reference script computes one gene per invocation, driven there by
run_variability_batch.sh + a manifest) and is not itself part of the ported
algorithm.

--taxid_tsv is a hard requirement of compute_variability.py (used for its
per-clade breakdown, <gene>.clade_entropy.tsv) — pass params.tax_id through
unchanged, whatever its schema. A 2-column tax_id file (this pipeline's own
auto-generated one, bin/generate_taxid_map.py) parses to zero clade rows
(load_taxonomy skips any line with fewer than 5 columns) rather than
erroring — CT_ACCUMULATION's own consumer never reads the clade file anyway
(subworkflows/CT_ACCUMULATION/local/src/aggregation/concatenate.py only
reads <gene>.entropy.tsv's `position`/`variability` columns) — but a richer
5-column taxid.tsv (tax_id, species, family, rank, name_class — the format
this pipeline's own fixtures already use) gets the full per-clade output.

compute_variability.py only understands FASTA (its own hand-rolled
read_fasta(), not Biopython) and derives the gene name via
`basename(prot_ali).replace('.fa', '')` — verbatim, including that it
replaces every '.fa' substring, not just a trailing extension. This driver
symlinks each alignment file to a temporary `<gene>.fa` before invoking it,
so gene-name extraction is correct regardless of the real alignment
filename/extension (--ali_format must be "fasta"; this path doesn't support
other formats since the reference script doesn't).

Usage
-----
    compute_alignment_entropy.py --alignment-dir <dir> --output-dir <dir> \
        --taxid-tsv <tax_id_file> [--family-order-tsv <tsv>] \
        [--alpha 1.0] [--beta 1.0] [--gamma 1.0]
"""

import argparse
import os
import subprocess
import sys
import tempfile

BIN_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--taxid-tsv", required=True,
                         help="Required by compute_variability.py itself (per-clade breakdown)")
    parser.add_argument("--family-order-tsv", default=None)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=1.0)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    files = sorted(
        f for f in os.listdir(args.alignment_dir)
        if os.path.isfile(os.path.join(args.alignment_dir, f))
    )
    if not files:
        print(f"Error: no alignment files found in {args.alignment_dir}", file=sys.stderr)
        sys.exit(1)

    n_ok = 0
    with tempfile.TemporaryDirectory() as tmp_dir:
        for fname in files:
            gene = os.path.splitext(fname)[0]
            src = os.path.abspath(os.path.join(args.alignment_dir, fname))
            link = os.path.join(tmp_dir, f"{gene}.fa")
            if not os.path.exists(link):
                os.symlink(src, link)

            cmd = [
                sys.executable, os.path.join(BIN_DIR, "compute_variability.py"),
                "--prot_ali", link,
                "--taxid_tsv", args.taxid_tsv,
                "--out_dir", args.output_dir,
                "--var_subdir", "",
                "--alpha", str(args.alpha),
                "--beta", str(args.beta),
                "--gamma", str(args.gamma),
            ]
            if args.family_order_tsv:
                cmd += ["--family_order_tsv", args.family_order_tsv]

            result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"WARN: compute_variability.py failed for {gene}: {result.stderr}",
                      file=sys.stderr)
                continue
            n_ok += 1

    print(f"Wrote {n_ok}/{len(files)} entropy files to {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
