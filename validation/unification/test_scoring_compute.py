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


def _stage(tmp_path):
    for f in ("scoring_compute.R", "aa_grouping.R"):
        shutil.copy(SRC / f, tmp_path / f)
    subprocess.run([sys.executable, str(SRC / "observed_core_scores.py"), "--input", str(GOLD / "filtered_discovery.tsv"),
                    "--positions-out", str(tmp_path / "core_pos.tsv"), "--genes-out", str(tmp_path / "core_genes.tsv")],
                   check=True, capture_output=True)


def _r(tmp_path, *extra, core=True):
    args = ["Rscript", "scoring_compute.R", "--postproc", str(GOLD / "filtered_discovery.tsv")]
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
                          "caas_sum": (row.CAAS_score - 1e-15) * 2, "n_schemes": 2} for i in range(3)])
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
