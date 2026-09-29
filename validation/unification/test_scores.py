"""core.scores: position mean, side collapse and size_adj_max, checked against the scores the R pipeline wrote."""
import math
import random
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.scores import (  # noqa: E402
    collapse_sides, direction_values, gene_scores, position_score, position_sum, size_adj_max,
)


def test_position_score_does_not_depend_on_fill_order():
    a = {"GS1": 0.1, "US": 0.3, "GS3": 0.7}
    b = dict(reversed(list(a.items())))
    assert position_sum(a) == position_sum(b)
    assert position_score(a) == pytest.approx((0.3 + 0.7 + 0.1) / 3)


def test_position_sum_is_correctly_rounded():
    """A naive left-to-right sum of these gives 0.0; the exact sum is 1.0."""
    assert position_sum({"US": 1e16, "GS4": 1.0, "GS3": -1e16}) == (1.0, 3)


def test_position_score_skips_missing_and_is_none_without_values():
    assert position_score({"US": 0.5, "GS4": float("nan"), "GS3": None}) == 0.5
    assert position_score({"US": None}) is None
    assert position_score({}) is None


def test_collapse_sides_counts_a_both_position_once():
    assert collapse_sides({"top": 0.4, "bottom": 0.9}) == {"all": 0.9, "top": 0.4, "bottom": 0.9}
    assert collapse_sides({"none": 0.2}) == {"all": 0.2}
    assert collapse_sides({"top": None}) == {}


def test_size_adj_max_closed_form_and_ties():
    pool = [0.1, 0.2, 0.2, 0.5, 0.9]
    assert size_adj_max([0.2, 0.05], pool) == pytest.approx((3 / 5) ** 2)   # ties at the max count as <=
    assert size_adj_max([0.95], pool) == 1.0
    assert size_adj_max([0.01], pool) == 0.0


def test_size_adj_max_is_none_when_there_is_nothing_to_score():
    assert size_adj_max([], [0.1, 0.2]) is None
    assert size_adj_max([float("nan"), None], [0.1, 0.2]) is None
    assert size_adj_max([0.3], []) is None


def test_gene_without_a_direction_has_no_score_for_it():
    pools = {"all": [0.1, 0.5, 0.9], "top": [0.1, 0.9], "bottom": [0.5]}
    got = gene_scores([{"bottom": 0.5}, {"bottom": 0.4}], pools)
    assert got["top"] is None
    assert got["bottom"] == pytest.approx(1.0 ** 2)
    assert got["all"] == pytest.approx((2 / 3) ** 2)


def test_size_adj_max_is_monotone_in_the_maximum():
    rng = random.Random(2)
    pool = sorted(rng.random() for _ in range(200))
    for n in (1, 3, 8):
        lo = size_adj_max([0.3] * n, pool)
        hi = size_adj_max([0.3] * (n - 1) + [0.6], pool)
        assert hi >= lo


# ── against the scores written by scoring_compute.R ──────────────────────────

def _r_reference(scoring_dir):
    pos = pd.read_csv(scoring_dir / "position_scores.tsv", sep="\t", keep_default_na=False, na_values=["NA", ""])
    genes = pd.read_csv(scoring_dir / "gene_scores.tsv", sep="\t", keep_default_na=False, na_values=["NA", ""])
    pos = pos.dropna(subset=["CAAS_score"])
    by_pos = {}
    for r in pos.itertuples():
        by_pos.setdefault((r.Gene, r.Position), {})[r.side] = r.CAAS_score
    pools = {d: sorted(v) for d, v in direction_values(by_pos.values()).items()}
    per_gene = {}
    for (g, _p), sides in by_pos.items():
        per_gene.setdefault(g, []).append(sides)
    return per_gene, pools, genes


@pytest.mark.parametrize("scoring_dir", [
    HERE / "golden/pepc_c4_complete",
    HERE / "cancer_b0_toy/neoplasia_prevalence_toy_complete/scoring",
], ids=["pepc", "toy"])
def test_gene_scores_match_the_r_pipeline(scoring_dir):
    if not (scoring_dir / "gene_scores.tsv").exists():
        pytest.skip(f"{scoring_dir} not available")
    per_gene, pools, genes = _r_reference(scoring_dir)
    checked = 0
    for r in genes.itertuples():
        got = gene_scores(per_gene[r.Gene], pools)
        for col, d in (("gene_caas_score", "all"), ("gene_caas_score_top", "top"), ("gene_caas_score_bottom", "bottom")):
            want = getattr(r, col)
            if pd.isna(want):
                assert got[d] is None, (r.Gene, col)
            else:
                assert got[d] == pytest.approx(want, abs=1e-12), (r.Gene, col)
                checked += 1
    assert checked > 0


# ── the null's gene x cycle scores ───────────────────────────────────────────

def test_null_gene_cycle_scores_use_core_scores_and_na_for_empty_cells(tmp_path):
    import csv
    import gzip

    from src.utils.gene_wrapper import _finalize_perm_scores

    d = tmp_path / "detail"
    d.mkdir()
    header = ["Gene", "cycle", "Position", "caap_group", "asr_path_score", "n_detected", "clust", "side"]
    rows = {  # (gene, position, side, [scheme scores]) in cycle c1; c2 has nothing
        "gA": [(1, "bottom", [0.2, 0.4]), (2, "bottom", [0.8])],
        "gB": [(1, "bottom", [0.6]), (5, "top", [0.3, 0.5, 0.7])],
    }
    schemes = ["US", "GS4", "GS3"]
    for g, ps in rows.items():
        with gzip.open(d / f"{g}.tsv.gz", "wt", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            for pos, side, vals in ps:
                for grp, v in zip(schemes, vals):
                    w.writerow([g, "c1", pos, grp, v, 1, 0, side])
    (tmp_path / "out").mkdir()
    _finalize_perm_scores(d, tmp_path / "out", ["c1", "c2"])
    got = pd.read_csv(tmp_path / "out/gene_cycle_scores.tsv", sep="\t", keep_default_na=True)
    got = got.set_index(["Gene", "cycle"])

    pos_scores = {("gA", 1): {"bottom": 0.3}, ("gA", 2): {"bottom": 0.8},
                  ("gB", 1): {"bottom": 0.6}, ("gB", 5): {"top": 0.5}}
    pools = {k: sorted(v) for k, v in direction_values(pos_scores.values()).items()}
    per_gene = {"gA": [pos_scores[("gA", 1)], pos_scores[("gA", 2)]],
                "gB": [pos_scores[("gB", 1)], pos_scores[("gB", 5)]]}
    for g, plist in per_gene.items():
        want = gene_scores(plist, pools)
        for d_, col in (("all", "global_caas"), ("top", "top_caas"), ("bottom", "bottom_caas")):
            v = got.loc[(g, "c1"), col]
            assert (pd.isna(v) and want[d_] is None) or v == pytest.approx(want[d_], abs=1e-12), (g, col)
    assert pd.isna(got.loc[("gA", "c1"), "top_caas"])           # gA has no top position
    assert got.loc[("gA", "c2"), ["global_caas", "top_caas", "bottom_caas"]].isna().all()   # empty cycle
    assert (got.loc[("gA", "c2"), ["global_asr", "top_asr", "bottom_asr"]] == 0).all()


# ── ties at the maximum ──────────────────────────────────────────────────────

def test_size_adj_max_counts_values_within_rounding_noise_of_the_maximum_as_ties():
    """Means that are equal in exact arithmetic can differ by an ulp; they must count as ties."""
    m = 0.3
    pool = sorted([0.1, m - 3e-17, m, math.nextafter(m, 1.0), math.nextafter(math.nextafter(m, 1.0), 1.0), 0.9])
    assert size_adj_max([m], pool) == pytest.approx(5 / 6)   # the three near-ties and 0.1 are <= max, 0.9 is not


def test_size_adj_max_keeps_real_differences_apart():
    pool = [0.1, 0.3, 0.3 + 1e-9, 0.9]
    assert size_adj_max([0.3], pool) == pytest.approx(2 / 4)


def test_size_adj_max_is_invariant_to_ulp_noise_in_the_inputs():
    rng = random.Random(4)
    base = [round(rng.random(), 3) for _ in range(300)]
    noisy = [math.nextafter(v, 1.0) if rng.random() < .5 else v for v in base]
    for m in base[:40]:
        assert size_adj_max([m], sorted(base)) == size_adj_max([m], sorted(noisy))
