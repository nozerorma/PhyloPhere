#!/usr/bin/env python3
"""Per-side permulation-null aggregation (_build_cycle_score_pools /
_finalize_perm_scores / _finalize_perm_pos_pval).

A per-side shard (a "both" position = two rows, one core_s each) is scored per
direction, and the global pool takes ONE max-deduped entry per position
(T3-doc §12), never both side rows.

Run: python -m pytest test_perm_side_split.py
"""
import bisect
import csv
import gzip
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # .../local
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import src.utils.gene_wrapper as gw  # noqa: E402

LEGACY = ["Gene", "cycle", "Position", "caap_group",
          "asr_path_score", "n_detected", "ct", "cb", "clust"]
SIDED = LEGACY + ["side"]


def _shard(dir_path: Path, gene: str, rows, fields):
    with gzip.open(dir_path / f"{gene}.tsv.gz", "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(dict(zip(fields, r)))


def _hist(rows, nd_idx=5, cyc_idx=1, pos_idx=2, grp_idx=3):
    h = {}
    seen = set()
    for r in rows:
        k = (r[cyc_idx], r[pos_idx], r[grp_idx])
        if k in seen:
            continue
        seen.add(k)
        h.setdefault(r[cyc_idx], {}).setdefault(r[nd_idx], 0)
        h[r[cyc_idx]][r[nd_idx]] += 1
    return h


def test_sided_shard_max_dedup():
    # G1 pos 10 is "both" -> two rows (top core 0.30, bottom core 0.70).
    # G1 pos 11 is top-only (0.20). One cycle, one scheme.
    rows = [("G1", "c1", 10, "US", 0.30, 2, 1, 0, 0, "top"),
            ("G1", "c1", 10, "US", 0.70, 2, 0, 1, 0, "bottom"),
            ("G1", "c1", 11, "US", 0.20, 1, 1, 0, 0, "top")]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "perm_pos_detail"
        d.mkdir()
        _shard(d, "G1", rows, SIDED)
        rl = gw.build_percent_rank_lookup(_hist(rows))
        pools = gw._build_cycle_score_pools(d, rl)
        # top pool: pos10 top 0.30, pos11 top 0.20
        assert sorted(pools["c1"]["top"]) == [0.20, 0.30]
        # bottom pool: pos10 bottom 0.70 only
        assert sorted(pools["c1"]["bottom"]) == [0.70]
        # global: ONE entry per position = its best side -> max(0.30,0.70)=0.70, and 0.20
        assert sorted(pools["c1"]["all"]) == [0.20, 0.70]

        gw._finalize_perm_scores(d, Path(td), cycle_tags=["c1"], rank_lookup=rl)
        scored = {}
        with open(Path(td) / "gene_cycle_scores.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                scored[(row["Gene"], row["cycle"])] = row
        g = scored[("G1", "c1")]
        # global_caas: gene G1 has 2 positions, max sides {0.70, 0.20}; F over the
        # cycle's `all` pool {0.20, 0.70} at max 0.70 -> 2/2 ; ^2 = 1.0
        assert abs(float(g["global_caas"]) - 1.0) < 1e-12, g
        # top_caas: G1 top scores {0.30, 0.20}; pool `top` {0.20,0.30}; max 0.30 ->
        # 2/2 ; ^2 = 1.0
        assert abs(float(g["top_caas"]) - 1.0) < 1e-12, g
        # bottom_caas: G1 bottom {0.70}; pool `bottom` {0.70}; 1/1 ^1 = 1.0
        assert abs(float(g["bottom_caas"]) - 1.0) < 1e-12, g


def test_sided_perm_pos_pval_two_rows():
    rows = [("G1", "c1", 10, "US", 0.30, 3, 1, 0, 0, "top"),
            ("G1", "c1", 10, "US", 0.70, 3, 0, 1, 0, "bottom"),
            ("G1", "c2", 10, "US", 0.30, 3, 1, 0, 0, "top"),
            ("G1", "c2", 10, "US", 0.70, 3, 0, 1, 0, "bottom"),
            ("G1", "c3", 10, "US", 0.30, 3, 1, 0, 0, "top"),
            ("G1", "c3", 10, "US", 0.70, 3, 0, 1, 0, "bottom"),
            ("G1", "c1", 11, "US", 0.20, 1, 1, 0, 0, "top")]
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "perm_pos_detail"
        d.mkdir()
        _shard(d, "G1", rows, SIDED)
        gw._finalize_perm_pos_pval(d, Path(td), cycle_tags=["c1", "c2", "c3"])
        got = {}
        with open(Path(td) / "perm_pos_pval.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                got[(row["Position"], row["side"])] = row
        # pos 10 "both" -> two rows, same n_detected / pos_perm_p
        assert ("10", "top") in got and ("10", "bottom") in got
        assert got[("10", "top")]["n_detected"] == got[("10", "bottom")]["n_detected"] == "3"
        assert got[("10", "top")]["pos_perm_p"] == got[("10", "bottom")]["pos_perm_p"]
        assert ("11", "top") in got and ("11", "bottom") not in got


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    ok = True
    for fn in fns:
        try:
            fn()
            print("PASS ", fn.__name__)
        except AssertionError as e:
            ok = False
            print("FAIL ", fn.__name__, "--", e)
    print("\n" + ("ALL TESTS PASSED" if ok else "SOME TESTS FAILED"))
    sys.exit(0 if ok else 1)
