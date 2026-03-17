#!/usr/bin/env python3
"""
collect_gene_sets.py
────────────────────
Collect significant gene sets from CT accumulation and CT postprocessing
outputs, and write consolidated gene lists for TOP and BOTTOM directions.

Sources integrated
──────────────────
TOP genes (≥1 source required):
  1. Accumulation CSV (convergent_caap_top):
     accumulation_convergent_caap_top_aggregated_results.csv
     → ActualCount_convergent_caap_top > 0
     AND PValueEmpirical_convergent_caap_top_FDR < fdr_threshold
  2. Postproc TXT: special_union_us_nondiv_and_us_gs_cases_change_side_top_significant.txt
     → gene names listed directly (no filtering needed, already significant)

BOTTOM genes (same logic):
  3. Accumulation CSV (convergent_caap_bottom):
     accumulation_convergent_caap_bottom_aggregated_results.csv
  4. Postproc TXT: special_union_us_nondiv_and_us_gs_cases_change_side_bottom_significant.txt

Usage
-----
    python collect_gene_sets.py \\
        [--accumulation_top  accumulated_top.csv] \\
        [--accumulation_bottom accumulated_bottom.csv] \\
        [--postproc_top      postproc_top.txt] \\
        [--postproc_bottom   postproc_bottom.txt] \\
        [--fdr               0.05] \\
        --out_top    gene_set_top.txt \\
        --out_bottom gene_set_bottom.txt

All source arguments are optional; at least one source per direction should be
supplied or the corresponding output file will be empty (a warning is issued).
"""

import argparse
import csv
import os
import sys


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_accumulation_csv(path: str, cat: str, fdr: float) -> set:
    """
    Parse a CT accumulation aggregated_results CSV and return genes above threshold.

    cat  : category string embedded in column names, e.g. 'convergent_caap_top'
    """
    count_col = f"ActualCount_{cat}"
    fdr_col   = f"PValueEmpirical_{cat}_FDR"

    genes = set()
    if not path or not os.path.isfile(path):
        return genes

    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            # Normalise column names (strip whitespace)
            fieldnames = [f.strip() for f in (reader.fieldnames or [])]
            reader.fieldnames = fieldnames

            if count_col not in fieldnames or fdr_col not in fieldnames:
                sys.stderr.write(
                    f"WARNING collect_gene_sets: expected columns [{count_col}, {fdr_col}] "
                    f"not found in '{path}'. Available: {fieldnames}\n"
                )
                return genes

            for row in reader:
                try:
                    count = float(row.get(count_col, 0) or 0)
                    fdr_v = float(row.get(fdr_col,   1) or 1)
                except ValueError:
                    continue
                if count > 0 and fdr_v < fdr:
                    gene_id = (row.get("GeneID") or row.get("Gene", "")).strip()
                    if gene_id:
                        genes.add(gene_id)
    except Exception as exc:
        sys.stderr.write(f"WARNING collect_gene_sets: could not read '{path}': {exc}\n")

    return genes


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
    parser.add_argument("--accumulation_top",    default="",
                        help="accumulation_convergent_caap_top_aggregated_results.csv")
    parser.add_argument("--accumulation_bottom", default="",
                        help="accumulation_convergent_caap_bottom_aggregated_results.csv")
    parser.add_argument("--postproc_top",        default="",
                        help="special_union_..._top_significant.txt")
    parser.add_argument("--postproc_bottom",     default="",
                        help="special_union_..._bottom_significant.txt")
    parser.add_argument("--fdr",                 type=float, default=0.05,
                        help="FDR threshold for accumulation filtering (default 0.05)")
    parser.add_argument("--out_top",    required=True, help="Output: gene_set_top.txt")
    parser.add_argument("--out_bottom", required=True, help="Output: gene_set_bottom.txt")
    args = parser.parse_args()

    # ── TOP ──────────────────────────────────────────────────────────────────
    top_genes: set = set()
    top_genes |= read_accumulation_csv(
        args.accumulation_top, "convergent_caap_top", args.fdr
    )
    top_genes |= read_txt_gene_list(args.postproc_top)

    if not top_genes:
        sys.stderr.write(
            "WARNING collect_gene_sets: no TOP genes collected. "
            "Check --accumulation_top and --postproc_top paths and FDR threshold.\n"
        )

    # ── BOTTOM ───────────────────────────────────────────────────────────────
    bottom_genes: set = set()
    bottom_genes |= read_accumulation_csv(
        args.accumulation_bottom, "convergent_caap_bottom", args.fdr
    )
    bottom_genes |= read_txt_gene_list(args.postproc_bottom)

    if not bottom_genes:
        sys.stderr.write(
            "WARNING collect_gene_sets: no BOTTOM genes collected. "
            "Check --accumulation_bottom and --postproc_bottom paths and FDR threshold.\n"
        )

    write_gene_list(args.out_top,    top_genes,    "TOP")
    write_gene_list(args.out_bottom, bottom_genes, "BOTTOM")

    # Summary
    sys.stderr.write(
        f"INFO collect_gene_sets: TOP={len(top_genes)}, BOTTOM={len(bottom_genes)} genes\n"
    )


if __name__ == "__main__":
    main()
