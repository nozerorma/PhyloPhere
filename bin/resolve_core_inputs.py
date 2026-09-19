#!/usr/bin/env python3
"""
resolve_core_inputs.py  —  Resolve tax_id / gene_ensembl_file before invoking
Nextflow, auto-generating whichever is left blank.

Why this runs outside main.nf: Nextflow (25.x) enforces single-assignment on
each params key — a params.tax_id = ... inside workflow{} that runs after
conf/common.config's own params.tax_id = params.tax_id ?: "" is silently
ignored ("`params.tax_id` is defined multiple times -- Assignments following
the first are ignored"), confirmed empirically, not assumed. So auto-
generation can't happen invisibly inside the pipeline for the ~15 places
that read params.tax_id/params.gene_ensembl_file directly — it has to
resolve to a real value BEFORE it reaches the `--tax_id`/`--gene_ensembl_file`
CLI flags. This script is that step; the GUI's generated run scripts call it
automatically, and a manual CLI run should too (see below).

Usage
-----
    resolve_core_inputs.py --outdir <outdir> \
        [--tree <tree.nwk>] [--alignment <alignment_dir>] \
        [--tax-id <existing_tax_id_path>] [--gene-ensembl-file <existing_path>]

Prints two lines to stdout, shell-sourceable:
    TAX_ID=<path or empty>
    GENE_ENSEMBL_FILE=<path or empty>

Example (manual CLI use)
-------------------------
    eval "$(python3 bin/resolve_core_inputs.py --outdir out --tree tree.nwk --alignment align/)"
    nextflow main.nf ... --tax_id "$TAX_ID" --gene_ensembl_file "$GENE_ENSEMBL_FILE"

A generation failure (missing ete3/pybiomart, network outage, no exact
matches) leaves the corresponding line empty rather than aborting — the
pipeline's own required-file checks then report it the same as if the user
had left the flag blank.
"""

import argparse
import os
import subprocess
import sys

BIN_DIR = os.path.dirname(os.path.abspath(__file__))


def resolve_tax_id(outdir: str, tree: str) -> str:
    core_dir = os.path.join(outdir, "core_inputs")
    os.makedirs(core_dir, exist_ok=True)
    output = os.path.join(core_dir, "tax_id_generated.tsv")
    unresolved = os.path.join(core_dir, "tax_id_unresolved.tsv")
    result = subprocess.run(
        [sys.executable, os.path.join(BIN_DIR, "generate_taxid_map.py"),
         "--tree", tree, "--output", output, "--unresolved", unresolved],
        stderr=subprocess.PIPE, text=True,
    )
    sys.stderr.write(result.stderr)
    if result.returncode != 0 or not os.path.exists(output):
        print("Warning: tax_id auto-generation failed; leaving --tax_id blank.",
              file=sys.stderr)
        return ""
    return output


def resolve_gene_ensembl_file(outdir: str, alignment_dir: str) -> str:
    core_dir = os.path.join(outdir, "core_inputs")
    os.makedirs(core_dir, exist_ok=True)
    genes = sorted({
        os.path.splitext(f)[0]
        for f in os.listdir(alignment_dir)
        if os.path.isfile(os.path.join(alignment_dir, f))
    })
    if not genes:
        print(f"Warning: no genes found in {alignment_dir}; leaving "
              "--gene_ensembl_file blank.", file=sys.stderr)
        return ""
    gene_list_file = os.path.join(core_dir, "gene_list.txt")
    with open(gene_list_file, "w") as fh:
        fh.write("\n".join(genes) + "\n")
    output = os.path.join(core_dir, "gene_ensembl_generated.tsv")
    unresolved = os.path.join(core_dir, "gene_ensembl_unresolved.txt")
    result = subprocess.run(
        [sys.executable, os.path.join(BIN_DIR, "generate_ensembl_mapping.py"),
         "--gene-list", gene_list_file, "--output", output, "--unresolved", unresolved],
        stderr=subprocess.PIPE, text=True,
    )
    sys.stderr.write(result.stderr)
    if result.returncode != 0 or not os.path.exists(output):
        print("Warning: gene_ensembl_file auto-generation failed; leaving "
              "--gene_ensembl_file blank.", file=sys.stderr)
        return ""
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tree", default="")
    parser.add_argument("--alignment", default="")
    parser.add_argument("--tax-id", default="", help="Existing --tax_id value, if any")
    parser.add_argument("--gene-ensembl-file", default="",
                         help="Existing --gene_ensembl_file value, if any")
    args = parser.parse_args()

    tax_id = args.tax_id
    if not tax_id and args.tree:
        tax_id = resolve_tax_id(args.outdir, args.tree)

    gene_ensembl_file = args.gene_ensembl_file
    if not gene_ensembl_file and args.alignment:
        gene_ensembl_file = resolve_gene_ensembl_file(args.outdir, args.alignment)

    print(f"TAX_ID={tax_id}")
    print(f"GENE_ENSEMBL_FILE={gene_ensembl_file}")


if __name__ == "__main__":
    main()
