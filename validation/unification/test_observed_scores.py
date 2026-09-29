"""observed_core_scores.py: the observed position and gene scores come from core.scores, and match the R pipeline's."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CLI = ROOT / "subworkflows/SCORING/local/src/observed_core_scores.py"
NA = dict(keep_default_na=False, na_values=["NA", ""])


def _run(tmp_path, inp):
    pos, genes = tmp_path / "pos.tsv", tmp_path / "genes.tsv"
    r = subprocess.run([sys.executable, str(CLI), "--input", str(inp), "--positions-out", str(pos),
                        "--genes-out", str(genes)], capture_output=True, text=True)
    return r, pos, genes


def _write(tmp_path, rows):
    p = tmp_path / "fd.tsv"
    pd.DataFrame(rows, columns=["Gene", "Position", "caap_group", "side", "asr_path_score"]).to_csv(
        p, sep="\t", index=False)
    return p


def test_scope_and_means(tmp_path):
    rows = [("gA", 1, "US", "top", 0.2), ("gA", 1, "GS4", "top", 0.4), ("gA", 1, "OTHER", "top", 0.99),
            ("gA", 2, "US", "bottom", 0.8), ("gA", 3, "US", "top", "NA"),
            ("gB", 1, "US", "bottom", 0.6)]
    r, pos, genes = _run(tmp_path, _write(tmp_path, rows))
    assert r.returncode == 0, r.stderr[-500:]
    p = pd.read_csv(pos, sep="\t", **NA).set_index(["Gene", "Position", "side"])
    assert p.loc[("gA", 1, "top"), "CAAS_score"] == pytest.approx(0.3)      # OTHER scheme is out of scope
    assert pd.isna(p.loc[("gA", 3, "top"), "CAAS_score"])                   # no numeric score, kept as NA
    g = pd.read_csv(genes, sep="\t", **NA).set_index("Gene")
    # pool all = {0.3, 0.8, 0.6}; gA max 0.8 with n=2 positions scored
    assert g.loc["gA", "gene_caas_score"] == pytest.approx(1.0 ** 2)
    assert g.loc["gB", "gene_caas_score"] == pytest.approx((2 / 3) ** 1)
    assert pd.isna(g.loc["gB", "gene_caas_score_top_all"])                  # gB has no top position


def test_a_scheme_repeated_within_a_position_is_an_error(tmp_path):
    rows = [("gA", 1, "US", "top", 0.2), ("gA", 1, "US", "top", 0.4)]
    r, _, _ = _run(tmp_path, _write(tmp_path, rows))
    assert r.returncode != 0 and "duplicate" in r.stderr.lower()


@pytest.mark.parametrize("fx", [
    (HERE / "golden/pepc_c4_complete/filtered_discovery.tsv", HERE / "golden/pepc_c4_complete"),
    (HERE / "cancer_b0_toy/neoplasia_prevalence_toy_complete/postproc/gene_filtering/filtered_discovery.tsv",
     HERE / "cancer_b0_toy/neoplasia_prevalence_toy_complete/scoring"),
], ids=["pepc", "toy"])
def test_matches_the_scores_written_by_r(tmp_path, fx):
    inp, ref = fx
    if not inp.exists() or not (ref / "position_scores.tsv").exists():
        pytest.skip("fixture not available")
    r, pos, genes = _run(tmp_path, inp)
    assert r.returncode == 0, r.stderr[-500:]
    got = pd.read_csv(pos, sep="\t", **NA, float_precision="round_trip").set_index(["Gene", "Position", "side"])
    want = pd.read_csv(ref / "position_scores.tsv", sep="\t", **NA).set_index(["Gene", "Position", "side"])
    assert set(got.index) == set(want.index)
    d = (got["CAAS_score"] - want["CAAS_score"].reindex(got.index)).abs()
    assert d.max() < 1e-12

    g = pd.read_csv(genes, sep="\t", **NA).set_index("Gene")
    wg = pd.read_csv(ref / "gene_scores.tsv", sep="\t", **NA).set_index("Gene")
    for got_col, want_col in (("gene_caas_score", "gene_caas_score"),
                              ("gene_caas_score_top_all", "gene_caas_score_top_all"),
                              ("gene_caas_score_bottom_all", "gene_caas_score_bottom_all")):
        a, b = g[got_col].reindex(wg.index), wg[want_col]
        assert (a.isna() == b.isna()).all(), got_col
        ok = a.notna()
        # R compared ties exactly; the core counts values within TIE_TOL as ties, so a gene may differ
        # where rounding noise moved R's count. Everything else agrees.
        diff = (a[ok] - b[ok]).abs()
        assert (diff > 1e-12).mean() < 0.02, (got_col, int((diff > 1e-12).sum()), int(ok.sum()))
