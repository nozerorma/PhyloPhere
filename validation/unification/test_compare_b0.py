"""Plumbing test for compare_b0.py on a synthetic two-gene run (observed == b_0 by construction)."""
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HARNESS = Path(__file__).with_name("compare_b0.py")


def _w(df, path, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if gz else None)


@pytest.fixture
def run(tmp_path):
    r = tmp_path / "run"
    # positions: G1@10 top (survivor), G1@11 top (cluster-flagged), G2@5 bottom (survivor)
    disc = pd.DataFrame({
        "gene": ["G1", "G1", "G2"], "mode": "CAAP", "caap_group": "US",
        "trait": ["traitfile_H1.tab"] * 3, "position": [10, 11, 5],
        "caas": ["AB/CC", "AB/CC", "DE/FF"], "amino_encoded": ["AB/CC", "AB/CC", "DE/FF"],
    })
    _w(disc, r / "caastools/discovery.tab")
    pd_disc = disc.assign(cycle="b_0~H1")[["cycle", "gene", "caap_group", "position", "caas", "amino_encoded"]]
    _w(pd_disc, r / "caas_permulation/perm_disc/G.perm_replay.discovery.output")

    master = pd.DataFrame({
        "gene": ["G1", "G1", "G2"], "msa_pos": [10, 11, 5], "caap_group": "US",
        "side": ["top", "top", "bottom"], "asr_path_score": [0.5, 0.25, 0.75]})
    (r / "ct_disambiguation").mkdir(parents=True)
    master.to_csv(r / "ct_disambiguation/caas_convergence_master.csv", index=False)
    _w(pd.DataFrame({"Gene": ["G1", "G2"], "Position": [10, 5], "caap_group": "US", "side": ["top", "bottom"]}),
       r / "postproc/gene_filtering/filtered_discovery.tsv")
    _w(pd.DataFrame({"Gene": ["G1", "G2"], "Position": [10, 5], "side": ["top", "bottom"], "CAAS_score": [0.5, 0.75]}),
       r / "scoring/position_scores.tsv")
    _w(pd.DataFrame({"Gene": ["G1", "G2"], "gene_caas_score": [0.4, 0.9],
                     "gene_caas_score_top_all": [0.4, None], "gene_caas_score_bottom_all": [None, 0.9]}),
       r / "scoring/gene_scores.tsv")

    b0 = r / "caas_permulation/b0"
    for g, rows in {"G1": [(10, 0.5, 0), (11, 0.25, 1)], "G2": [(5, 0.75, 0)]}.items():
        side = "top" if g == "G1" else "bottom"
        _w(pd.DataFrame([{"Gene": g, "cycle": "b_0", "Position": p, "caap_group": "US", "asr_path_score": s,
                          "n_detected": 1, "clust": c, "side": side} for p, s, c in rows]),
           b0 / f"perm_pos_detail/{g}.tsv.gz", gz=True)
    _w(pd.DataFrame({"Gene": ["G1", "G2"], "Position": [10, 5], "side": ["top", "bottom"], "cycle": "b_0",
                     "caas_score": [0.5, 0.75], "n_schemes": 1}), b0 / "perm_pos_cycle_caas.tsv.gz", gz=True)
    _w(pd.DataFrame({"Gene": ["G1", "G2"], "cycle": "b_0", "global_asr": 0.0, "top_asr": 0.0, "bottom_asr": 0.0,
                     "global_caas": [0.4, 0.9], "top_caas": [0.4, None], "bottom_caas": [None, 0.9]}),
       b0 / "gene_cycle_scores.tsv")
    return r


def _run(r):
    p = subprocess.run([sys.executable, str(HARNESS), "--run", str(r)], capture_output=True, text=True)
    return p.returncode, p.stdout


def test_identical_passes(run):
    rc, out = _run(run)
    assert rc == 0, out
    assert out.count("PASS") == 5


def test_value_perturbation_fails_only_c(run):
    m = pd.read_csv(run / "ct_disambiguation/caas_convergence_master.csv")
    m.loc[0, "asr_path_score"] += 1e-9  # above the 1e-12 tolerance
    m.to_csv(run / "ct_disambiguation/caas_convergence_master.csv", index=False)
    rc, out = _run(run)
    assert rc == 1
    assert "[C] FAIL" in out and "[A] PASS" in out and "[D] PASS" in out


def test_missing_discovery_row_fails_a(run):
    p = run / "caas_permulation/perm_disc/G.perm_replay.discovery.output"
    d = pd.read_csv(p, sep="\t")
    d.iloc[:-1].to_csv(p, sep="\t", index=False)
    rc, out = _run(run)
    assert rc == 1 and "[A] FAIL" in out


def test_removed_unit_changes_survivors_b(run):
    # a removed (b_0, US, G2) unit drops G2 from the b_0 survivors, so B must fail
    _w(pd.DataFrame({"cycle": ["b_0"], "caap_group": ["US"], "Gene": ["G2"]}),
       run / "caas_permulation/b0/removed_units.tsv")
    rc, out = _run(run)
    assert rc == 1 and "[B] FAIL" in out


def test_a_zero_where_the_observed_gene_score_is_na_fails_e(run):
    g = pd.read_csv(run / "caas_permulation/b0/gene_cycle_scores.tsv", sep="\t")
    g["top_caas"] = g["top_caas"].fillna(0.0)          # the null's earlier convention for an empty direction
    g.to_csv(run / "caas_permulation/b0/gene_cycle_scores.tsv", sep="\t", index=False)
    rc, out = _run(run)
    assert rc == 1 and "[E] FAIL" in out and out.count("PASS") == 4


def test_c_reports_bitwise_differences(run):
    m = pd.read_csv(run / "ct_disambiguation/caas_convergence_master.csv")
    m.loc[0, "asr_path_score"] = 0.5 + 1e-16          # within tolerance, not the same bits
    m.to_csv(run / "ct_disambiguation/caas_convergence_master.csv", index=False)
    rc, out = _run(run)
    assert rc == 0 and '"n_bitwise_different": 1' in out
