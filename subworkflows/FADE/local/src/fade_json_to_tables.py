#!/usr/bin/env python3
# fade_json_to_tables.py — Reduce one HyPhy FADE JSON to the gene-level row and the significant-site rows of its gene.
# PhyloPhere | subworkflows/FADE/local/src/

"""
Reads the Bayes factor of every (site, target amino acid) of one ``<gene>.<direction>.FADE.json`` and writes two small tab-delimited
tables, so that the JSON itself does not have to be kept:

    --genes-out   one row: direction, gene, n_sites_tested, n_sig_site_aa, max_bf, top_target_aa, top_site, sigcount_<AA> (20 columns)
    --sites-out   one row per site with at least one BF >= threshold: direction, gene, site (1-based), position (0-based), max_bf, target_aa

Layout of the JSON as read here: ``MLE.headers`` (list of [column, description]) and ``MLE.content[AA][partition]`` (one list of
numbers per site, in the order of the headers). Only the ``BayesFactor`` column is used.

* n_sites_tested   the largest number of sites of any partition
* n_sig_site_aa    number of (site, AA) pairs with BF >= threshold
* max_bf           the largest finite BF, with its site (top_site, 1-based) and amino acid (top_target_aa); both NA below the threshold
* sigcount_<AA>    number of significant sites of each target amino acid
* sites table      the largest BF of each significant site over the 20 amino acids, and the amino acid that gives it
                   (the first one in alphabetical order when two tie)

A JSON without MLE content, or without any finite BF, gives no row (exit status 3): the gene is skipped.

Imported by: subworkflows/FADE/local/src/run_hyphy_fade_batch.sh, subworkflows/FADE/fade_run.nf (FADE_RUN)
Usage: fade_json_to_tables.py --json F --gene G --direction D --bf-threshold T --genes-out P --sites-out P
       fade_json_to_tables.py --print-header genes|sites
"""

import argparse
import json
import math
import sys

AA_LETTERS = "ACDEFGHIKLMNPQRSTVWY"
GENES_COLUMNS = (["direction", "gene", "n_sites_tested", "n_sig_site_aa", "max_bf", "top_target_aa", "top_site"]
                 + [f"sigcount_{aa}" for aa in AA_LETTERS])
SITES_COLUMNS = ["direction", "gene", "site", "position", "max_bf", "target_aa"]


def _bf_column(headers):
    """0-based index of the Bayes factor column (the fourth one when no header names it)."""
    for i, h in enumerate(headers or []):
        if "bayesfactor" in str(h[0]).lower():
            return i
    return 3


def _num(x):
    return float(x) if isinstance(x, (int, float)) and not isinstance(x, bool) else None


def reduce_json(js, gene, direction, threshold):
    """(gene row, site rows) of one parsed FADE JSON, or None when it holds no usable result."""
    mle = js.get("MLE") or {}
    content = mle.get("content")
    if not content:
        return None
    bf_col = _bf_column(mle.get("headers"))
    n_sites_tested, n_sig = 0, 0
    best_bf, best_site, best_aa = -math.inf, None, None
    sig_counts = {aa: 0 for aa in AA_LETTERS}
    sig_site_bf, sig_site_aa = {}, {}          # 1-based site -> largest BF over the amino acids and the amino acid that gives it
    for aa in AA_LETTERS:                      # alphabetical: the first amino acid wins a tie
        by_partition = content.get(aa)
        if not by_partition:
            continue
        for rows in by_partition.values():
            if not rows:
                continue
            n_sites_tested = max(n_sites_tested, len(rows))
            bf = [_num(r[bf_col]) if len(r) > bf_col else None for r in rows]
            for site, v in enumerate(bf, 1):
                if v is not None and v >= threshold:
                    n_sig += 1
                    sig_counts[aa] += 1
                    if site not in sig_site_bf or v > sig_site_bf[site]:
                        sig_site_bf[site], sig_site_aa[site] = v, aa
            finite = [(v, site) for site, v in enumerate(bf, 1) if v is not None and math.isfinite(v)]
            if finite:
                v, site = max(finite, key=lambda t: (t[0], -t[1]))     # the first site of the largest BF
                if v > best_bf:
                    best_bf, best_site, best_aa = v, site, aa
    if not math.isfinite(best_bf):
        return None
    hit = best_bf >= threshold
    row = [direction, gene, n_sites_tested, n_sig, best_bf, best_aa if hit else "NA", best_site if hit else "NA"]
    row += [sig_counts[aa] for aa in AA_LETTERS]
    sites = [[direction, gene, s, s - 1, sig_site_bf[s], sig_site_aa[s]] for s in sig_site_bf]
    return row, sites


def _line(values):
    """A tab-delimited line; floats carry 15 significant digits, as R writes them."""
    return "\t".join(format(v, ".15g") if isinstance(v, float) else str(v) for v in values) + "\n"


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--print-header":      # the header of a table without rows
        sys.stdout.write(_line({"genes": GENES_COLUMNS, "sites": SITES_COLUMNS}[sys.argv[2]]))
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--json", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--direction", required=True)
    ap.add_argument("--bf-threshold", type=float, default=100.0)
    ap.add_argument("--genes-out", required=True)
    ap.add_argument("--sites-out", required=True)
    a = ap.parse_args()
    try:
        with open(a.json) as fh:
            js = json.load(fh)
        reduced = reduce_json(js, a.gene, a.direction, a.bf_threshold)
    except (OSError, ValueError) as e:
        print(f"[fade_json_to_tables] {a.gene} ({a.direction}): {e}", file=sys.stderr)
        return 3
    if reduced is None:
        print(f"[fade_json_to_tables] {a.gene} ({a.direction}): no MLE content with a finite Bayes factor", file=sys.stderr)
        return 3
    row, sites = reduced
    with open(a.genes_out, "w") as fh:
        fh.write(_line(GENES_COLUMNS))
        fh.write(_line(row))
    with open(a.sites_out, "w") as fh:
        fh.write(_line(SITES_COLUMNS))
        for s in sites:
            fh.write(_line(s))
    return 0


if __name__ == "__main__":
    sys.exit(main())
