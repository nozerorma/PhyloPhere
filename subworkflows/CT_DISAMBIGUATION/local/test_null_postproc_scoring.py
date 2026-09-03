#!/usr/bin/env python3
"""Gap B integration test: CT_POSTPROC filtering of the CAAS permulation null.

Builds a synthetic perm_pos_detail.tsv.gz, runs pass B0 (cycle-aware gene
filter) + pass B (scoring), and checks:
  * the dubious gene-unit (count outlier + clustered) is removed for its cycle;
  * the removed unit contributes global_caas == 0 for that cycle;
  * a kept gene's global_caas equals size_adj_max recomputed over the FILTERED
    per-cycle pool (i.e. the null CAAS_score uses the same algebra as
    scoring_compute.R section 2f/2g/4a, on the post-filter pool).

Run: python test_null_postproc_scoring.py
"""
import bisect
import csv
import gzip
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "src"))

import src.utils.gene_wrapper as gw

DETAIL_FIELDS = ["Gene", "cycle", "Position", "caap_group",
                 "asr_path_score", "n_detected", "ct", "cb", "clust"]


def _rows():
    # Pass A writes the detail file one gene at a time, so it is gene-contiguous;
    # emit it the same way here (all of a gene's cycles together).
    r = []
    c1n = {"A": 2, "B": 2, "C": 3, "D": 2}
    c2n = {"A": 2, "B": 3, "C": 2, "D": 2}
    for g in ("A", "B", "C", "D"):
        for i in range(c1n[g]):
            r.append((g, "c1", 100 + i, "US", 0.4 + 0.02 * i, c1n[g], 1, 0, 0))
        for i in range(c2n[g]):
            r.append((g, "c2", 100 + i, "US", 0.5 + 0.02 * i, c2n[g], 0, 1, 0))
    for i in range(60):                       # E: c1 count outlier, all clustered
        r.append(("E", "c1", 500 + i, "US", 0.9, 60, 1, 0, 1))
    return r


def _write_detail(path, rows):
    with gzip.open(path, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=DETAIL_FIELDS, delimiter="\t")
        w.writeheader()
        for g, c, p, grp, asr, nd, ct, cb, cl in rows:
            w.writerow(dict(zip(DETAIL_FIELDS, (g, c, p, grp, asr, nd, ct, cb, cl))))


def main():
    rows = _rows()
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        detail = td / "perm_pos_detail.tsv.gz"
        _write_detail(detail, rows)

        hist = {}
        for _g, c, _p, _grp, _asr, nd, *_ in rows:
            hist.setdefault(c, {}).setdefault(nd, 0)
            hist[c][nd] += 1
        rank_lookup = gw.build_percent_rank_lookup(hist)

        removed = gw._cycle_gene_removal_from_detail(
            detail, gene_lengths={}, mode="dubious",
            iqr_multiplier=3.0, extreme_percentile=0.99)
        assert ("c1", "US", "E") in removed, removed
        assert not any(g == "E" and c == "c2" for c, _grp, g in removed), removed

        gw._finalize_perm_scores(detail, td, cycle_tags=["c1", "c2"],
                                 rank_lookup=rank_lookup, removed=removed)

        scored = {}
        with open(td / "gene_cycle_scores.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                scored[(row["Gene"], row["cycle"])] = row
        assert float(scored[("E", "c1")]["global_caas"]) == 0.0, scored[("E", "c1")]

        # Recompute the FILTERED c1/all pool and size_adj_max for gene C by hand.
        pool = []
        for g, c, p, grp, asr, nd, ct, cb, cl in rows:
            if c != "c1" or (c, grp, g) in removed:
                continue
            phen = 1.0 - rank_lookup["c1"].get(nd, 0.0)
            pool.append(phen * asr)           # one scheme per position -> mean == value
        pool.sort()
        c_vals = sorted(
            (1.0 - rank_lookup["c1"].get(nd, 0.0)) * asr
            for g, c, p, grp, asr, nd, *_ in rows if c == "c1" and g == "C")
        expect = (bisect.bisect_right(pool, max(c_vals)) / len(pool)) ** len(c_vals)
        got = float(scored[("C", "c1")]["global_caas"])
        assert abs(got - expect) < 1e-12, (got, expect)
        print(f"PASS  E removed from c1; C c1 global_caas={got:.6f} == size_adj_max over filtered pool")
        print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
