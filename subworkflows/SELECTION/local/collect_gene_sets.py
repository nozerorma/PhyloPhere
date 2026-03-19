#!/usr/bin/env python3
"""
collect_gene_sets.py
────────────────────
Collect significant gene sets from CT postprocessing outputs and write
consolidated gene lists for TOP and BOTTOM directions.

Sources integrated
──────────────────
TOP genes:
  Postproc TXT: special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt
  → gene names listed directly (no filtering needed, already significant)

BOTTOM genes:
  Postproc TXT: special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt

Usage
-----
    python collect_gene_sets.py \\
        [--postproc_top      postproc_top.txt] \\
        [--postproc_bottom   postproc_bottom.txt] \\
        --out_top    gene_set_top.txt \\
        --out_bottom gene_set_bottom.txt

Postproc source arguments are optional; missing directions produce empty
output files with a warning.
"""

import argparse
import os
import sys


def read_txt_gene_list(path: str) -> set:
    """Read a plain-text file with one gene name per line."""
    genes = set()
    if not path or not os.path.isfile(path):
        return genes
    try:
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith("#"):
                    genes.add(g)
    except Exception as exc:
        sys.stderr.write(f"WARNING collect_gene_sets: could not read '{path}': {exc}\n")
    return genes


def write_gene_list(path: str, genes: set, label: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for g in sorted(genes):
            fh.write(g + "\n")
    sys.stderr.write(
        f"INFO collect_gene_sets: wrote {len(genes)} {label} genes to '{path}'\n"
    )


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Collect TOP / BOTTOM significant gene sets for selection analyses"
    )
    parser.add_argument("--postproc_top",        default="",
                        help="special_union_..._top_significant.txt")
    parser.add_argument("--postproc_bottom",     default="",
                        help="special_union_..._bottom_significant.txt")
    parser.add_argument("--out_top",    required=True, help="Output: gene_set_top.txt")
    parser.add_argument("--out_bottom", required=True, help="Output: gene_set_bottom.txt")
    args = parser.parse_args()

    # ── TOP ──────────────────────────────────────────────────────────────────
    top_genes: set = set()
    top_genes |= read_txt_gene_list(args.postproc_top)

    if not top_genes:
        sys.stderr.write(
            "WARNING collect_gene_sets: no TOP genes collected. "
            "Check --postproc_top path.\n"
        )

    # ── BOTTOM ───────────────────────────────────────────────────────────────
    bottom_genes: set = set()
    bottom_genes |= read_txt_gene_list(args.postproc_bottom)

    if not bottom_genes:
        sys.stderr.write(
            "WARNING collect_gene_sets: no BOTTOM genes collected. "
            "Check --postproc_bottom path.\n"
        )

    write_gene_list(args.out_top,    top_genes,    "TOP")
    write_gene_list(args.out_bottom, bottom_genes, "BOTTOM")

    # Summary
    sys.stderr.write(
        f"INFO collect_gene_sets: TOP={len(top_genes)}, BOTTOM={len(bottom_genes)} genes\n"
    )


if __name__ == "__main__":
    main()
