"""scoring_compute.R integrates the core scores: it reads them, and compares against the null with a tie tolerance."""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "subworkflows/SCORING/local/src"
GOLD = HERE / "golden/pepc_c4_complete"
NA = dict(keep_default_na=False, na_values=["NA", ""])

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _stage(tmp_path, positions=None):
    """Scripts and core tables of the frozen PEPC run; `positions` keeps only the first N positions of its postproc table."""
    for f in ("scoring_compute.R", "aa_grouping.R"):
        shutil.copy(SRC / f, tmp_path / f)
    source = GOLD / "filtered_discovery.tsv"
    if positions:
        df = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        keep = df[["Gene", "Position"]].drop_duplicates().head(positions)
        source = tmp_path / "small_discovery.tsv"
        df.merge(keep).to_csv(source, sep="\t", index=False)
    (tmp_path / "postproc.txt").write_text(str(source))
    subprocess.run([sys.executable, str(SRC / "observed_core_scores.py"), "--input", str(source),
                    "--positions-out", str(tmp_path / "core_pos.tsv"), "--genes-out", str(tmp_path / "core_genes.tsv")],
                   check=True, capture_output=True)


def _r(tmp_path, *extra, core=True):
    args = ["Rscript", "scoring_compute.R", "--postproc", (tmp_path / "postproc.txt").read_text()]
    if core:
        args += ["--core_positions", "core_pos.tsv", "--core_genes", "core_genes.tsv"]
    return subprocess.run(args + list(extra), cwd=tmp_path, capture_output=True, text=True)


def test_scores_are_read_from_the_core_and_match_the_r_pipeline(tmp_path):
    _stage(tmp_path)
    r = _r(tmp_path)
    assert r.returncode == 0, r.stderr[-800:]
    got = pd.read_csv(tmp_path / "position_scores.tsv", sep="\t", **NA).set_index(["Gene", "Position", "side"])
    want = pd.read_csv(GOLD / "position_scores.tsv", sep="\t", **NA).set_index(["Gene", "Position", "side"])
    assert set(got.index) == set(want.index)
    assert (got["CAAS_score"] - want["CAAS_score"].reindex(got.index)).abs().max() < 1e-12
    g = pd.read_csv(tmp_path / "gene_scores.tsv", sep="\t", **NA).set_index("Gene")
    wg = pd.read_csv(GOLD / "gene_scores.tsv", sep="\t", **NA).set_index("Gene")
    for c in ("gene_caas_score", "gene_caas_score_top", "gene_caas_score_bottom"):
        assert ((g[c] - wg[c].reindex(g.index)).abs().fillna(0) < 1e-12).all(), c
        assert (g[c].isna() == wg[c].reindex(g.index).isna()).all(), c


def test_missing_core_tables_stop_the_script(tmp_path):
    _stage(tmp_path)
    r = _r(tmp_path, core=False)
    assert r.returncode != 0 and "core_positions" in r.stderr


def test_null_values_within_rounding_noise_of_the_observed_score_count_as_ties(tmp_path):
    _stage(tmp_path)
    core = pd.read_csv(tmp_path / "core_pos.tsv", sep="\t", **NA)
    row = core[core.CAAS_score > 0].iloc[0]
    null = pd.DataFrame([{"Gene": row.Gene, "Position": row.Position, "side": row.side, "cycle": f"c{i}",
                          "caas_score": row.CAAS_score - 1e-15, "n_schemes": 2} for i in range(3)])
    null.to_csv(tmp_path / "null.tsv.gz", sep="\t", index=False)
    r = _r(tmp_path, "--caas_pos_cycle_caas", str(tmp_path / "null.tsv.gz"))
    assert r.returncode == 0, r.stderr[-800:]
    ps = pd.read_csv(tmp_path / "position_scores.tsv", sep="\t", **NA)
    p = ps[(ps.Gene == row.Gene) & (ps.Position == row.Position)]["p.emp"].iloc[0]
    assert p == pytest.approx(1.0)      # all 3 cycles tie or exceed: (3 + 1) / (3 + 1)


def test_tie_tolerance_is_the_same_constant_in_r_and_python():
    import re
    sys.path.insert(0, str(HERE.parents[1] / "subworkflows/CT_DISAMBIGUATION/local"))
    from src.core.scores import TIE_TOL
    m = re.search(r"^TIE_TOL <- ([0-9.eE+-]+)", (SRC / "scoring_compute.R").read_text(), re.M)
    assert m and float(m.group(1)) == TIE_TOL


def test_a_null_without_caas_score_is_rejected(tmp_path):
    _stage(tmp_path)
    core = pd.read_csv(tmp_path / "core_pos.tsv", sep="\t", **NA)
    row = core[core.CAAS_score > 0].iloc[0]
    pd.DataFrame([{"Gene": row.Gene, "Position": row.Position, "side": row.side, "cycle": "c1",
                   "caas_sum": 0.5, "n_schemes": 1}]).to_csv(tmp_path / "old.tsv.gz", sep="\t", index=False)
    r = _r(tmp_path, "--caas_pos_cycle_caas", str(tmp_path / "old.tsv.gz"))
    assert r.returncode != 0 and "caas_score" in r.stderr


def _empty_null(tmp_path):
    """The null tables of N = 0, written by the producers themselves: header-only tables and the RDS built from them."""
    (tmp_path / "empty_detail").mkdir()
    subprocess.run([sys.executable, "./reaggregate_perm_scores.py", "--detail", str(tmp_path / "empty_detail"), "--output-dir",
                    str(tmp_path / "null"), "--empty-null"], cwd=SRC.parents[2] / "CT_DISAMBIGUATION/local", check=True, capture_output=True)
    subprocess.run(["Rscript", str(SRC / "scoring_caas_perms.R"), "--gene-cycle-scores", str(tmp_path / "null/gene_cycle_scores.tsv"),
                    "--output", str(tmp_path / "null/caas_perms.rds")], check=True, capture_output=True)
    return tmp_path / "null"


@pytest.mark.parametrize("positions", [None, 4], ids=["many_positions", "few_positions"])
def test_an_empty_null_leaves_every_null_based_value_undefined_and_keeps_the_observed_scores(tmp_path, positions):
    _stage(tmp_path, positions)
    assert _r(tmp_path).returncode == 0
    without = pd.read_csv(tmp_path / "position_scores.tsv", sep="\t", **NA)
    gene_without = pd.read_csv(tmp_path / "gene_scores.tsv", sep="\t", **NA)
    null = _empty_null(tmp_path)
    r = _r(tmp_path, "--caas_pos_cycle_caas", str(null / "perm_pos_cycle_caas.tsv.gz"), "--caas_perms", str(null / "caas_perms.rds"))
    assert r.returncode == 0, r.stderr[-800:]
    ps = pd.read_csv(tmp_path / "position_scores.tsv", sep="\t", **NA)
    assert ps[["p.emp", "p.adj_bh", "p.adj_sam"]].isna().all().all()
    gs = pd.read_csv(tmp_path / "gene_scores.tsv", sep="\t", **NA)
    assert gs[[c for c in gs.columns if "pperm" in c]].isna().all().all()
    # the observed side does not depend on the null
    assert ps["CAAS_score"].equals(without["CAAS_score"]) and gs["gene_caas_score"].equals(gene_without["gene_caas_score"])
    assert "p.emp" in (r.stdout + r.stderr) and "no null" in (r.stdout + r.stderr).lower()
    assert "no --caas_pos_cycle_caas provided" not in r.stdout      # the table was provided; it is the null that is empty
