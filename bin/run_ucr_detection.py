#!/usr/bin/env python3
# run_ucr_detection.py — Batch driver running detect_ucr.py on every per-gene entropy table of a directory.
# PhyloPhere | bin/

"""
RunUcrDetection: loops detect_ucr.py (one gene per call) over the <gene>.entropy.tsv files of
a directory; the detection logic is entirely in detect_ucr.py.

Called by:  RUN_UCR_DETECTION Nextflow process (subworkflows/ENRICHMENT/ucr_generation.nf → run_ucr_detection.py)
Inputs:     --entropy-dir   directory of <gene>.entropy.tsv (compute_alignment_entropy.py);
                            <gene>.clade_entropy.tsv files are ignored
            any other argument is forwarded to detect_ucr.py for every gene
            (--abs_threshold, --rel_zscore, --window_size ...; see its --help)
Outputs:    <output-dir>/<gene>.ucr.tsv for the genes with at least one region (detect_ucr.py writes
            no file otherwise); a count of the genes processed on stderr. A failure of one gene
            is reported on stderr and the loop continues.

Usage:
    run_ucr_detection.py --entropy-dir <dir> --output-dir <dir> [detect_ucr.py options...]
"""

import argparse
import glob
import os
import subprocess
import sys

BIN_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--entropy-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args, passthrough = parser.parse_known_args()

    os.makedirs(args.output_dir, exist_ok=True)

    entropy_files = sorted(glob.glob(os.path.join(args.entropy_dir, "*.entropy.tsv")))
    entropy_files = [f for f in entropy_files if not f.endswith(".clade_entropy.tsv")]
    if not entropy_files:
        print(f"Error: no *.entropy.tsv files found in {args.entropy_dir}", file=sys.stderr)
        sys.exit(1)

    n_ok = 0
    for path in entropy_files:
        gene = os.path.basename(path).replace(".entropy.tsv", "")
        out_tsv = os.path.join(args.output_dir, f"{gene}.ucr.tsv")
        cmd = [
            sys.executable, os.path.join(BIN_DIR, "detect_ucr.py"),
            "--entropy_tsv", path, "--gene", gene, "--out_tsv", out_tsv,
        ] + passthrough
        result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"WARN: detect_ucr.py failed for {gene}: {result.stderr}", file=sys.stderr)
            continue
        n_ok += 1

    print(f"Ran UCR detection for {n_ok}/{len(entropy_files)} genes -> {args.output_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
