"""The observed gene filter and the null pass B run the same post-processing at the same grain."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.postproc import GeneUnit, gene_removal  # noqa: E402
from src.utils.gene_wrapper import (  # noqa: E402
    _build_cycle_score_pools, _cycle_gene_removal_from_detail,
)

FILTER_GENES = ROOT / "subworkflows/CT_POSTPROC/local/src/filter_caas_genes.py"


def _table(tmp_path):
    """Nine quiet genes of 4 positions (the IQR threshold is 4, so g1 is not an outlier despite its 3-position
    train) plus gX (30 contiguous positions, all in a train)."""
    rows, disc = [], []
    for i in range(1, 10):
        pos = [1000 * i + j for j in (0, 20, 40, 60)] if i != 1 else [100, 101, 102, 300]
        rows += [(f"g{i}", p, "US") for p in pos]
    rows += [("gX", p, "US") for p in range(1, 31)]
    disc = [("gX", p, "US") for p in range(1, 31)] + [("g1", p, "US") for p in (100, 101, 102)]
    pd.DataFrame(rows, columns=["Gene", "Position", "caap_group"]).to_csv(tmp_path / "disc.tsv", sep="\t", index=False)
    cl = pd.DataFrame(disc, columns=["Gene", "Position", "caap_group"])
    cl["clustering_flag"] = "Discarded"
    cl.to_csv(tmp_path / "clusters.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": [f"g{i}" for i in range(1, 10)] + ["gX"], "length": 1000}).to_csv(
        tmp_path / "len.tsv", sep="\t", index=False)
    units = [GeneUnit("b_0", "US", g, n, g in {"g1", "gX"})
             for g, n in pd.DataFrame(rows).groupby(0)[1].nunique().items()]
    return units


def _run(tmp_path, *extra):
    cmd = [sys.executable, str(FILTER_GENES), "-i", str(tmp_path / "disc.tsv"), "-l", str(tmp_path / "len.tsv"),
           "-c", str(tmp_path / "clusters.tsv"), "-m", "dubious", "-o", str(tmp_path / "out.tsv"),
           "-s", str(tmp_path / "summary.tsv"), "-g", str(tmp_path / "stats.tsv"), *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_observed_gene_filter_removes_what_the_core_removes_on_the_pooled_rows(tmp_path):
    units = _table(tmp_path)
    r = _run(tmp_path)
    assert r.returncode == 0, r.stderr[-600:]
    expected = gene_removal(units, {u.gene: 1000.0 for u in units}, "dubious")
    assert set(expected) == {("b_0", "US", "gX")}
    summary = pd.read_csv(tmp_path / "summary.tsv", sep="\t")
    assert set(zip(summary["caap_group"], summary["Gene"], summary["category"])) == {("US", "gX", "Dubious")}
    out = pd.read_csv(tmp_path / "out.tsv", sep="\t")
    assert "gX" not in set(out["Gene"])
    # trains are kept unless --remove-clusters is given
    assert {100, 101, 102} <= set(out.loc[out["Gene"] == "g1", "Position"])
    stats = pd.read_csv(tmp_path / "stats.tsv", sep="\t")
    assert set(stats["hyp_id"]) == {"ALL"} and {"threshold_extreme", "threshold_dubious", "category"} <= set(stats.columns)
    assert stats.loc[stats["Gene"] == "gX", "category"].item() == "Dubious"


def test_observed_remove_clusters_drops_train_positions(tmp_path):
    _table(tmp_path)
    r = _run(tmp_path, "--remove-clusters")
    assert r.returncode == 0, r.stderr[-600:]
    out = pd.read_csv(tmp_path / "out.tsv", sep="\t")
    assert set(out.loc[out["Gene"] == "g1", "Position"]) == {300}


def _detail(tmp_path, clust):
    d = tmp_path / "detail"
    d.mkdir()
    with gzip.open(d / "g1.tsv.gz", "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", "cycle", "Position", "caap_group", "asr_path_score", "n_detected", "clust", "side"])
        for pos, c in zip(range(1, 6), clust):
            w.writerow(["g1", "c1", pos, "US", 0.5 + pos / 10, 1, c, "top"])
    return d


def test_null_pool_keeps_train_rows_when_clusters_are_not_removed(tmp_path):
    d = _detail(tmp_path, [1, 1, 0, 0, 0])
    assert len(_build_cycle_score_pools(d)["c1"]["all"]) == 3
    assert len(_build_cycle_score_pools(d, remove_clusters=False)["c1"]["all"]) == 5


def test_null_gene_removal_agrees_with_core_on_detail_rows(tmp_path):
    d = tmp_path / "detail"
    d.mkdir()
    genes = {f"g{i}": 2 for i in range(1, 10)} | {"gX": 30}
    with gzip.open(d / "all.tsv.gz", "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", "cycle", "Position", "caap_group", "asr_path_score", "n_detected", "clust", "side"])
        for g, n in genes.items():
            for p in range(1, n + 1):
                w.writerow([g, "c1", p, "US", 0.5, 1, int(g == "gX"), "top"])
    got = _cycle_gene_removal_from_detail(d, {g: 1000.0 for g in genes}, "dubious", 3.0, 0.99)
    assert got == {("c1", "US", "gX")}
