"""N of the null is the number of cycles replayed, not the number of cycles that left a row: the roster comes from the labelings."""
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from test_reaggregate_empty import LOCAL, R_SCRIPT, _shard

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "subworkflows/SCORING/local/src"
NA = dict(sep="\t", keep_default_na=False, na_values=["NA", ""])
needs_r = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _labelings(path, tags):
    path.write_text("".join(f"{t}\tsp1,sp2\tsp3,sp4\n" for t in tags))
    return path


def _run(detail, out, *extra):
    return subprocess.run([sys.executable, "./reaggregate_perm_scores.py", "--detail", str(detail), "--output-dir", str(out), *extra],
                          cwd=LOCAL, capture_output=True, text=True)


def _two_cycle_detail(d, tags=("b_1", "b_2")):
    _shard(d, "G1", tags[0])
    _shard(d, "G2", tags[1])
    return d


@pytest.mark.parametrize("tags", [
    ["b_0", "b_1", "b_2", "b_3"],                                                            # plain labelings
    ["b_0~H1", "b_0~H2", "b_1~H1", "b_1~H2", "b_2~H1", "b_3~H1", "b_3~H2"],                  # one row per hypothesis
])
def test_the_roster_lists_every_replayed_cycle_without_the_real_labeling(tmp_path, tags):
    _two_cycle_detail(tmp_path / "d")
    r = _run(tmp_path / "d", tmp_path / "out", "--cycles-from", str(_labelings(tmp_path / "lab.tab", tags)))
    assert r.returncode == 0, r.stderr[-800:]
    assert (tmp_path / "out/cycle_roster.txt").read_text().split() == ["b_1", "b_2", "b_3"]
    assert "1 of 3" in r.stderr and "b_3" in r.stderr                                          # the cycle without a row is reported


def test_the_roster_keeps_the_grain_of_the_detail_even_when_it_names_labelings(tmp_path):
    _two_cycle_detail(tmp_path / "d", ("b_1~H1", "b_2~H7"))          # hypotheses not pooled: a detail row carries the tag of its labeling
    r = _run(tmp_path / "d", tmp_path / "out", "--cycles-from", str(_labelings(tmp_path / "lab.tab", ["b_0~H1", "b_1~H1", "b_1~H7", "b_2~H7", "b_3~H1"])))
    assert r.returncode == 0, r.stderr[-800:]
    assert (tmp_path / "out/cycle_roster.txt").read_text().split() == ["b_1~H1", "b_1~H7", "b_2~H7", "b_3~H1"] and "2 of 4" in r.stderr


def test_a_cycle_with_rows_that_is_not_in_the_labelings_is_an_error(tmp_path):
    _two_cycle_detail(tmp_path / "d")
    r = _run(tmp_path / "d", tmp_path / "out", "--cycles-from", str(_labelings(tmp_path / "lab.tab", ["b_0", "b_1"])))
    assert r.returncode != 0 and "b_2" in r.stderr and "not in" in r.stderr


def test_without_labelings_nothing_changes(tmp_path):
    _two_cycle_detail(tmp_path / "d")
    assert _run(tmp_path / "d", tmp_path / "out").returncode == 0 and not (tmp_path / "out/cycle_roster.txt").exists()


def _gcs(path, cycles):
    pd.DataFrame([{"Gene": "G1", "cycle": c, "global_asr": 0.0, "top_asr": 0.0, "bottom_asr": 0.0, "global_caas": 0.5, "top_caas": 0.5,
                   "bottom_caas": 0.0} for c in cycles]).to_csv(path, sep="\t", index=False)


def _rds_cycles(rds):
    out = subprocess.run(["Rscript", "-e", f'x <- readRDS("{rds}"); cat(colnames(x$caas_corStat_byrank$global), sep=",")'], capture_output=True, text=True)
    return out.stdout.split(",")


@needs_r
def test_the_rds_has_a_column_for_every_cycle_of_the_roster(tmp_path):
    _gcs(tmp_path / "gcs.tsv", ["b_1", "b_2"])
    (tmp_path / "roster.txt").write_text("b_1\nb_2\nb_3\n")
    r = subprocess.run(["Rscript", str(R_SCRIPT), "--gene-cycle-scores", str(tmp_path / "gcs.tsv"), "--cycles", str(tmp_path / "roster.txt"),
                        "--output", str(tmp_path / "with.rds")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    assert _rds_cycles(tmp_path / "with.rds") == ["b_1", "b_2", "b_3"]
    subprocess.run(["Rscript", str(R_SCRIPT), "--gene-cycle-scores", str(tmp_path / "gcs.tsv"), "--output", str(tmp_path / "without.rds")], check=True, capture_output=True)
    assert _rds_cycles(tmp_path / "without.rds") == ["b_1", "b_2"]


@needs_r
def test_a_roster_that_misses_a_cycle_with_rows_stops_the_build(tmp_path):
    _gcs(tmp_path / "gcs.tsv", ["b_1", "b_2"])
    (tmp_path / "roster.txt").write_text("b_1\n")
    r = subprocess.run(["Rscript", str(R_SCRIPT), "--gene-cycle-scores", str(tmp_path / "gcs.tsv"), "--cycles", str(tmp_path / "roster.txt"),
                        "--output", str(tmp_path / "x.rds")], capture_output=True, text=True)
    assert r.returncode != 0 and "b_2" in r.stderr


@needs_r
def test_p_emp_counts_the_cycles_of_the_roster_not_only_the_cycles_with_rows(tmp_path):
    from test_scoring_compute import _stage, _r
    _stage(tmp_path)
    core = pd.read_csv(tmp_path / "core_pos.tsv", **NA)
    row = core[core.CAAS_score > 0].iloc[0]
    null = pd.DataFrame([{"Gene": row.Gene, "Position": row.Position, "side": row.side, "cycle": c, "caas_score": row.CAAS_score + 1.0,
                          "n_schemes": 2} for c in ("b_1", "b_2")])                 # every cycle with rows beats the observed
    null.to_csv(tmp_path / "null.tsv.gz", sep="\t", index=False)
    pd.DataFrame([{"Gene": row.Gene, "cycle": c, "global_asr": 0.0, "top_asr": 0.0, "bottom_asr": 0.0, "global_caas": 1.0, "top_caas": 1.0,
                   "bottom_caas": 0.0} for c in ("b_1", "b_2")]).to_csv(tmp_path / "gcs.tsv", sep="\t", index=False)
    pvalues = {}
    for name, roster in (("two", ["b_1", "b_2"]), ("four", ["b_1", "b_2", "b_3", "b_4"])):
        (tmp_path / f"{name}.txt").write_text("\n".join(roster) + "\n")
        subprocess.run(["Rscript", str(R_SCRIPT), "--gene-cycle-scores", str(tmp_path / "gcs.tsv"), "--cycles", str(tmp_path / f"{name}.txt"),
                        "--output", str(tmp_path / f"{name}.rds")], check=True, capture_output=True)
        r = _r(tmp_path, "--caas_pos_cycle_caas", str(tmp_path / "null.tsv.gz"), "--caas_perms", str(tmp_path / f"{name}.rds"))
        assert r.returncode == 0, r.stderr[-500:]
        ps = pd.read_csv(tmp_path / "position_scores.tsv", **NA)
        pvalues[name] = ps[(ps.Gene == row.Gene) & (ps.Position == row.Position)]["p.emp"].iloc[0]
    assert pvalues["two"] == pytest.approx(3 / 3) and pvalues["four"] == pytest.approx(3 / 5)   # (k + 1) / (N + 1), k = 2 either way
