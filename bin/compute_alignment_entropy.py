#!/usr/bin/env python3
# compute_alignment_entropy.py — Batch driver computing per-gene Valdar variability from an alignment directory.
# PhyloPhere | bin/

"""
ComputeAlignmentEntropy: runs compute_variability.py once per alignment file of a directory,
so that the per-gene <gene>.entropy.tsv files exist without a prior ortholog_characterizator run.

compute_variability.py (same directory) computes one gene per call; this driver only adds the
batching. All numbers come from that script.

Called by:  COMPUTE_ALIGNMENT_ENTROPY Nextflow process (subworkflows/CT_ACCUMULATION/ctacc_run.nf →
            compute_alignment_entropy.py); consumed by CT_ACCUMULATION (accumulation_entropy_dir)
            and by UCR_GENERATION
Inputs:     --alignment-dir    directory of FASTA alignments, one file per gene
            --taxid-tsv        required by compute_variability.py for its per-clade output
            --family-order-tsv optional family → order table; --alpha/--beta/--gamma Valdar exponents
Outputs:    <output-dir>/<gene>.entropy.tsv, <gene>.clade_entropy.tsv, <gene>.fa (variable columns
            only); a count of the genes done on stderr

Notes:
  * compute_variability.py reads only FASTA and derives the gene name with
    basename(path).replace('.fa', ''), which replaces every '.fa' substring. Each alignment is
    therefore symlinked as <gene>.fa in a temporary directory, with <gene> the file name
    without its last extension.
  * Clade rows need a taxid table with at least 5 columns (species in column 2, family in
    column 3, name class "scientific name" in column 5). A 2-column table such as the one of
    generate_taxid_map.py gives an empty clade file, with no error. The consumers
    (CT_ACCUMULATION, detect_ucr.py, aggregate_ucr.py) read only <gene>.entropy.tsv.

Usage:
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
