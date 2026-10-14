#!/usr/bin/env python3
# fade_merge_tables.py — Merge the per-batch FADE tables of one direction into the published tables.
# PhyloPhere | subworkflows/FADE/local/src/

"""
Concatenates the tables written by fade_json_to_tables.py (one pair per FADE task) and keeps the rows of one direction. Writes:

    fade_summary_<direction>.tsv      gene, n_sites_tested, n_sig_site_aa, max_bf, top_target_aa, top_site; one row per gene, by decreasing
                                      n_sig_site_aa and max_bf (read by SCORING and the gene lists)
    fade_gene_counts_<direction>.tsv  the same rows with the sigcount_<AA> columns (read by the report)
    fade_site_bf_<direction>.tsv      gene, position (0-based), max_bf, target_aa; one row per significant site, by decreasing max_bf
    fade_sites_<direction>.csv        the same rows, comma-delimited (read by SCORING and ENRICHMENT)

Every table has its header when no gene or site reached the file. Ties keep the order of the gene names, then the order of the sites
of the gene.

Imported by: subworkflows/FADE/fade_tables.nf (FADE_TABLES)
Usage: fade_merge_tables.py --direction D --genes F [F ...] --sites F [F ...]
"""

import argparse
import csv
import sys

SUMMARY = ["gene", "n_sites_tested", "n_sig_site_aa", "max_bf", "top_target_aa", "top_site"]
SITES = ["gene", "position", "max_bf", "target_aa"]


def _read(paths, direction):
    """Rows (dicts) of the tables of `paths` that belong to `direction`, in file order."""
    rows = []
    for p in paths:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["direction"] == direction:
                    rows.append(r)
    return rows


def _write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in header) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--direction", required=True)
    ap.add_argument("--genes", nargs="*", default=[])
    ap.add_argument("--sites", nargs="*", default=[])
    a = ap.parse_args()
    d = a.direction

    genes = sorted(_read(a.genes, d), key=lambda r: r["gene"])
    genes.sort(key=lambda r: (-int(r["n_sig_site_aa"]), -float(r["max_bf"])))       # stable: ties stay in gene order
    counts = [c for c in (genes[0] if genes else {}) if c.startswith("sigcount_")]
    if not counts:
        counts = [f"sigcount_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"]
    _write_tsv(f"fade_summary_{d}.tsv", SUMMARY, genes)
    _write_tsv(f"fade_gene_counts_{d}.tsv", SUMMARY + counts, genes)

    sites = sorted(_read(a.sites, d), key=lambda r: r["gene"])
    sites.sort(key=lambda r: -float(r["max_bf"]))
    _write_tsv(f"fade_site_bf_{d}.tsv", SITES, sites)
    with open(f"fade_sites_{d}.csv", "w") as fh:       # quoted strings and bare numbers, as write.csv writes them
        q = lambda x: '"' + str(x).replace('"', '""') + '"'
        fh.write(",".join(q(c) for c in SITES) + "\n")
        for r in sites:
            fh.write(f'{q(r["gene"])},{int(r["position"])},{r["max_bf"]},{q(r["target_aa"])}\n')
    print(f"[fade_merge_tables] {d}: {len(genes)} genes, {len(sites)} significant sites")
    return 0


if __name__ == "__main__":
    sys.exit(main())
