"""core.postproc: cluster trains and per-labeling gene removal, shared by the observed chain and the null."""
import itertools
import random
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.postproc import GeneUnit, ctrain, gene_removal, gene_unit_stats, train_flags  # noqa: E402


# ── trains ───────────────────────────────────────────────────────────────────

def _ctrain_bruteforce(positions, maxcaas, minlen):
    """Every index interval of the sorted unique positions, tested independently of the loop in ctrain."""
    u = sorted(set(positions))
    bad = set()
    for l, r in itertools.combinations_with_replacement(range(len(u)), 2):
        span = u[r] - u[l] + 1
        if span >= minlen and (r - l + 1) / span >= maxcaas:
            bad.update(u[l:r + 1])
    return sorted(bad)


def test_ctrain_hand_cases():
    assert ctrain([10, 11, 12, 50], 0.7, 3) == [10, 11, 12]
    assert ctrain([10, 30, 50], 0.7, 3) == []
    assert ctrain([10, 11], 0.7, 3) == []           # fewer positions than minlen
    assert ctrain([5, 5, 6, 7], 0.7, 3) == [5, 6, 7]  # duplicates collapse


def test_ctrain_matches_bruteforce_on_random_inputs():
    rng = random.Random(7)
    for _ in range(300):
        pos = [rng.randint(1, 40) for _ in range(rng.randint(0, 14))]
        maxcaas, minlen = rng.choice([0.5, 0.7, 0.9]), rng.choice([2, 3, 5])
        assert ctrain(pos, maxcaas, minlen) == _ctrain_bruteforce(pos, maxcaas, minlen)


def test_train_flags_subset_of_union():
    """Trains of a subset of the positions are trains of the union: adding positions only raises density."""
    rng = random.Random(11)
    for _ in range(200):
        union = {rng.randint(1, 50) for _ in range(rng.randint(3, 20))}
        sub = {p for p in union if rng.random() < 0.6}
        f_sub = train_flags({"k": sub}, 0.7, 3)["k"]
        f_union = train_flags({"k": union}, 0.7, 3)["k"]
        assert f_sub <= f_union


def test_train_flags_keeps_keys_independent():
    got = train_flags({("b_0", "US"): [1, 2, 3], ("b_0", "G2"): [1, 20, 40]}, 0.7, 3)
    assert got == {("b_0", "US"): {1, 2, 3}, ("b_0", "G2"): set()}


# ── gene removal ─────────────────────────────────────────────────────────────

def _units(labeling, group, counts, clustered=()):
    return [GeneUnit(labeling, group, g, n, g in clustered) for g, n in counts.items()]


LEN = {f"g{i}": 1000.0 for i in range(1, 12)}


def test_dubious_needs_iqr_outlier_and_a_train():
    counts = {f"g{i}": 2 for i in range(1, 10)} | {"g10": 40, "g11": 40}
    # g10 is an outlier with a train, g11 an outlier without one
    got = gene_removal(_units("b_0", "US", counts, clustered={"g10"}), LEN, "dubious")
    assert got == {("b_0", "US", "g10"): "Dubious"}


def test_dubious_is_strictly_above_the_threshold():
    """Equal counts give IQR 0: every unit sits ON the threshold, so none is an outlier even with trains."""
    counts = {f"g{i}": 7 for i in range(1, 8)}
    assert gene_removal(_units("b_0", "US", counts, clustered=set(counts)), LEN, "dubious") == {}


def test_extreme_is_top_density_and_ignores_genes_without_length():
    counts = {f"g{i}": 5 for i in range(1, 11)} | {"g11": 500}
    lens = dict(LEN)
    got = gene_removal(_units("b_0", "US", counts), lens, "extreme", extreme_percentile=0.9)
    assert set(got) == {("b_0", "US", "g11")} and set(got.values()) == {"Extreme"}
    # without an annotated length the gene has no density, so it is never removed
    del lens["g11"]
    assert gene_removal(_units("b_0", "US", counts), lens, "extreme", extreme_percentile=0.9) == {}


def test_dubious_pool_includes_genes_without_length():
    """Dubious counts positions only: an unannotated gene still calibrates the IQR pool."""
    counts = {f"g{i}": 2 for i in range(1, 6)} | {"gX": 30}
    lens = {g: 1000.0 for g in counts if g != "gX"}
    stats = {s["Gene"]: s for s in gene_unit_stats(_units("b_0", "US", counts, clustered={"gX"}), lens)}
    assert stats["gX"]["length"] is None and stats["gX"]["density"] is None
    assert stats["g1"]["threshold_dubious"] == pytest.approx(2.0)  # Q1 = Q3 = 2
    assert stats["gX"]["dubious"] is True                         # unannotated, still an outlier


def test_both_category_and_none_mode():
    counts = {f"g{i}": 2 for i in range(1, 10)} | {"g10": 60}
    u = _units("b_0", "US", counts, clustered={"g10"})
    assert gene_removal(u, LEN, "both", extreme_percentile=0.9) == {("b_0", "US", "g10"): "Both"}
    assert gene_removal(u, LEN, "none") == {}


def test_pools_are_per_labeling_and_group():
    """The same counts are an outlier in a quiet labeling and not in a noisy one."""
    quiet = {f"g{i}": 2 for i in range(1, 10)} | {"g10": 40}
    noisy = {f"g{i}": 2 * i for i in range(1, 10)} | {"g10": 40}
    u = _units("q", "US", quiet, {"g10"}) + _units("n", "US", noisy, {"g10"}) + _units("q", "G2", noisy, {"g10"})
    assert set(gene_removal(u, LEN, "dubious")) == {("q", "US", "g10")}


def test_thresholds_match_pandas_quantile():
    rng = random.Random(3)
    for _ in range(50):
        n = rng.randint(2, 40)
        counts = {f"g{i}": rng.randint(1, 60) for i in range(n)}
        lens = {g: rng.uniform(200, 5000) for g in counts}
        for k, q in [(1.5, 0.9), (3.0, 0.99)]:
            stats = gene_unit_stats(_units("b_0", "US", counts), lens, k, q)
            ser = pd.Series(counts, dtype=float)
            q1, q3 = ser.quantile(0.25), ser.quantile(0.75)
            dens = pd.Series({g: counts[g] / lens[g] * 100 for g in counts})
            assert stats[0]["threshold_dubious"] == pytest.approx(q3 + k * (q3 - q1))
            assert stats[0]["threshold_extreme"] == pytest.approx(dens.quantile(q))
            assert {s["Gene"] for s in stats if s["extreme"]} == set(dens[dens > dens.quantile(q)].index)
