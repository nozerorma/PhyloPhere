#!/usr/bin/env python3
# resolve_core_inputs.py — Resolve tax_id and gene_ensembl_file before Nextflow starts, generating the blank ones.
# PhyloPhere | bin/

"""
ResolveCoreInputs: generates the tax_id table (generate_taxid_map.py) and the gene Ensembl table
(generate_ensembl_mapping.py) that the user left blank, and prints the paths for the caller.

Why it runs outside main.nf: Nextflow keeps the first assignment of a params key and ignores
later ones, so a params.tax_id set inside the workflow after conf/common.config has assigned its
default would have no effect on the processes that read it. The value must exist before it
reaches the --tax_id / --gene_ensembl_file flags.

Called by:  the generated run script (gui/generation/templates/run_single.sh.j2), before nextflow;
            a manual run should do the same (see Usage)
Inputs:     --outdir              results directory; files go to <outdir>/core_inputs/
            --tree                species tree, for tax_id
            --alignment           alignment directory (gene = file name without extension), for gene_ensembl_file
            --tax-id, --gene-ensembl-file   values already given; they are returned unchanged
            --auto-generate-ensembl         required to generate gene_ensembl_file; also --ensembl-dataset, --ref-species
Outputs:    stdout, shell-sourceable:  TAX_ID=<path or empty>  and  GENE_ENSEMBL_FILE=<path or empty>
            files: tax_id_generated.tsv, tax_id_unresolved.tsv, gene_list.txt,
            gene_ensembl_generated.tsv, gene_ensembl_unresolved.txt

A failed generation (missing ete3 or requests, network outage, no exact match) leaves its
line empty instead of aborting; the pipeline's own checks of required files then report it as
for a blank flag.

Usage:
    resolve_core_inputs.py --outdir <outdir> \
        [--tree <tree.nwk>] [--alignment <alignment_dir>] [--auto-generate-ensembl] \
        [--tax-id <existing_tax_id_path>] [--gene-ensembl-file <existing_path>]

    eval "$(python3 bin/resolve_core_inputs.py --outdir out --tree tree.nwk --alignment align/)"
    nextflow main.nf ... --tax_id "$TAX_ID" --gene_ensembl_file "$GENE_ENSEMBL_FILE"
"""

import argparse
import os
import subprocess
import sys

BIN_DIR = os.path.dirname(os.path.abspath(__file__))


def resolve_tax_id(outdir: str, tree: str) -> str:
    """Run generate_taxid_map.py on the tree; returns the table path, or "" on failure."""
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


def resolve_gene_ensembl_file(outdir: str, alignment_dir: str,
                              dataset: str = "", ref_species: str = "") -> str:
    """Run generate_ensembl_mapping.py on the genes of the alignment directory; returns the table path, or "" on failure."""
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
    cmd = [
        sys.executable, os.path.join(BIN_DIR, "generate_ensembl_mapping.py"),
        "--gene-list", gene_list_file, "--output", output, "--unresolved", unresolved,
    ]
    if dataset:
        cmd.extend(["--dataset", dataset])
    if ref_species:
        cmd.extend(["--ref-species", ref_species])

    result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
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
    parser.add_argument("--auto-generate-ensembl", action="store_true", default=False,
                         help="Allow auto-generating gene_ensembl_file via BioMart if not provided")
    parser.add_argument("--ensembl-dataset", default="",
                         help="BioMart dataset name")
    parser.add_argument("--ref-species", default="Homo_sapiens",
                         help="Reference species name")
    args = parser.parse_args()

    tax_id = args.tax_id
    if not tax_id and args.tree:
        tax_id = resolve_tax_id(args.outdir, args.tree)

    gene_ensembl_file = args.gene_ensembl_file
    if not gene_ensembl_file and args.alignment and args.auto_generate_ensembl:
        gene_ensembl_file = resolve_gene_ensembl_file(
            args.outdir, args.alignment,
            dataset=args.ensembl_dataset, ref_species=args.ref_species
        )

    print(f"TAX_ID={tax_id}")
    print(f"GENE_ENSEMBL_FILE={gene_ensembl_file}")


if __name__ == "__main__":
    main()
