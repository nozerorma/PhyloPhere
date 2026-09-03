#!/usr/bin/env python3
"""Tests for null_postproc.py — cycle-aware CT_POSTPROC filtering of the CAAS null.

Run: python -m pytest test_null_postproc.py   (or: python test_null_postproc.py)
"""
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "null_postproc", Path(__file__).with_name("null_postproc.py"))
_np = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_np)
ctrain = _np.ctrain
clustering_discards = _np.clustering_discards
cycle_gene_removal = _np.cycle_gene_removal


def test_ctrain_matches_docstring_example():
    # filter_caas_clusters-param.py docstring: [10,11,12,50], maxcaas=0.7, minlen=3
    assert ctrain([10, 11, 12, 50], maxcaas=0.7, minlen=3) == [10, 11, 12]


def test_ctrain_below_minlen_is_noop():
    assert ctrain([5, 40], maxcaas=0.7, minlen=3) == []
    assert ctrain([1, 2], maxcaas=0.5, minlen=3) == []


def test_ctrain_sparse_positions_kept():
    assert ctrain([1, 20, 45, 90], maxcaas=0.7, minlen=3) == []


def test_ctrain_dense_block():
    # 5 consecutive positions: span 5, count 5, density 1.0
    assert ctrain([100, 101, 102, 103, 104], maxcaas=0.7, minlen=3) == [100, 101, 102, 103, 104]


def test_clustering_discards_per_group():
    out = clustering_discards(
        {"US": [10, 11, 12, 50], "GS4": [10, 60, 90]}, maxcaas=0.7, minlen=3
    )
    assert out["US"] == {10, 11, 12}
    assert out["GS4"] == set()


def test_cycle_gene_removal_dubious_needs_cluster_and_outlier():
    # cycle c1 / group US: E is a count outlier AND clustered -> dubious.
    # cycle c2 has the same count distribution but its outlier F is NOT
    # clustered -> kept (this is the "needs cluster" half of the assertion).
    rows = [
        ("c1", "US", "A", 2, False), ("c1", "US", "B", 2, False),
        ("c1", "US", "C", 3, False), ("c1", "US", "D", 2, False),
        ("c1", "US", "E", 200, True),
        ("c2", "US", "A", 2, False), ("c2", "US", "B", 2, False),
        ("c2", "US", "C", 3, False), ("c2", "US", "D", 2, False),
        ("c2", "US", "F", 200, False),
    ]
    rem = cycle_gene_removal(rows, gene_lengths={}, mode="dubious", iqr_multiplier=3.0)
    assert ("c1", "US", "E") in rem
    assert ("c2", "US", "F") not in rem


def test_cycle_gene_removal_is_per_cycle():
    # Same gene E: count-outlier+clustered in c1, ordinary in c2 -> only c1.
    rows = [
        ("c1", "US", "A", 2, False), ("c1", "US", "B", 2, False),
        ("c1", "US", "C", 3, False), ("c1", "US", "D", 2, False),
        ("c1", "US", "E", 200, True),
        ("c2", "US", "A", 20, False), ("c2", "US", "B", 22, False),
        ("c2", "US", "C", 19, False), ("c2", "US", "D", 20, False),
        ("c2", "US", "E", 21, True),
    ]
    rem = cycle_gene_removal(rows, gene_lengths={}, mode="dubious")
    assert ("c1", "US", "E") in rem
    assert ("c2", "US", "E") not in rem


def test_cycle_gene_removal_extreme_uses_density():
    rows = [
        ("c1", "US", "A", 5, False), ("c1", "US", "B", 5, False),
        ("c1", "US", "C", 5, False), ("c1", "US", "D", 5, False),
        ("c1", "US", "E", 6, False),
    ]
    # E is short -> very high density -> extreme; others long.
    lengths = {"A": 1000, "B": 1000, "C": 1000, "D": 1000, "E": 50}
    rem = cycle_gene_removal(rows, gene_lengths=lengths, mode="extreme",
                             extreme_percentile=0.99)
    assert ("c1", "US", "E") in rem
    assert ("c1", "US", "A") not in rem


def test_cycle_gene_removal_none_mode():
    rows = [("c1", "US", "E", 40, True), ("c1", "US", "A", 1, False)]
    assert cycle_gene_removal(rows, {}, mode="none") == set()


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    sys.exit(1 if failed else 0)
