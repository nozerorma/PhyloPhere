"""The length-stratified gene removal (gene_filter_mode = stratified): the fences of core.postproc.stratified_fences, the
decisions of gene_removal and gene_unit_stats, and filter_caas_genes.py on the observed table."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PIPELINE = Path(__file__).resolve().parents[1]
LOCAL = PIPELINE / "subworkflows/CT_DISAMBIGUATION/local"
FILTER = PIPELINE / "subworkflows/CT_POSTPROC/local/src/filter_caas_genes.py"
sys.path.insert(0, str(LOCAL))

from src.core.postproc import GeneUnit, gene_removal, gene_unit_stats, stratified_fences  # noqa: E402


def test_the_fences_are_computed_inside_the_length_strata():
    genes = ["a", "b", "c", "d"]
    counts = [1, 2, 3, 50]
    lengths = {"a": 100, "b": 200, "c": 300, "d": 400}
    far_near = stratified_fences(genes, counts, lengths, strata=2, k_far=3.0, k_near=1.5)
    # pool IQR = 14.75 - 1.75 = 13; stratum {1, 2}: Q3 = 1.75, IQR 0.5 -> scale 13; stratum {3, 50}: Q3 = 38.25, IQR 23.5
    assert far_near[0] == far_near[1] == pytest.approx((1.75 + 3 * 13, 1.75 + 1.5 * 13))
    assert far_near[2] == far_near[3] == pytest.approx((38.25 + 3 * 23.5, 38.25 + 1.5 * 23.5))


def test_a_gene_without_a_length_has_no_fence_and_is_never_removed():
    genes = ["a", "b", "nolen"]
    assert stratified_fences(genes, [1, 2, 900], {"a": 100, "b": 200}, strata=2)[2] == (None, None)
    units = [GeneUnit("b_0", g, n, True) for g, n in zip(genes, [1, 2, 900])]
    assert ("b_0", "nolen") not in gene_removal(units, {"a": 100, "b": 200}, mode="stratified")


def _pool(n_genes=60):
    """Counts that grow with the length, plus one far outlier, one near outlier with a train and one near outlier without."""
    lengths = {f"g{i:02d}": 100 + 20 * i for i in range(n_genes)}
    counts = {g: 1 + (i % 4) + i // 6 for i, g in enumerate(lengths)}          # 1 .. 13, rising with the length
    return lengths, counts


def test_categories_follow_the_two_fences_and_the_train_flag():
    lengths, counts = _pool()
    counts.update({"g05": 60, "g06": 60, "g30": 25})                          # g05 far + train, g06 far without train, g30 near
    trains = {"g05", "g30"}
    units = [GeneUnit("b_0", g, n, g in trains) for g, n in counts.items()]
    removed = gene_removal(units, lengths, mode="stratified")
    assert removed[("b_0", "g05")] == "Both"                                  # above the far fence, with a train
    assert removed[("b_0", "g06")] == "Extreme"                               # above the far fence, no train needed
    assert removed.get(("b_0", "g30")) in ("Dubious", "Both")                 # above the near fence, with a train
    assert not any(g in ("g00", "g01", "g59") for (_lab, g) in removed)       # ordinary genes stay


def test_the_stats_and_the_removal_agree_and_the_fences_are_per_gene():
    lengths, counts = _pool()
    counts["g06"] = 60
    units = [GeneUnit("b_0", g, n, g == "g05") for g, n in counts.items()]
    removed = gene_removal(units, lengths, mode="stratified")
    stats = gene_unit_stats(units, lengths, mode="stratified")
    flagged = {s["Gene"]: ("Both" if s["extreme"] and s["dubious"] else "Dubious" if s["dubious"] else "Extreme")
               for s in stats if s["extreme"] or s["dubious"]}
    assert flagged == {g: cat for (_lab, g), cat in removed.items()}
    assert len({s["threshold_extreme"] for s in stats if s["threshold_extreme"] is not None}) > 1   # one fence per stratum


def test_the_removal_does_not_depend_on_the_length_of_a_gene_by_itself():
    """A pool whose counts scale with the length: the long genes are not removed for being long."""
    lengths = {f"g{i:03d}": 100 + 10 * i for i in range(200)}
    counts = {g: 1 + (100 + 10 * i) // 80 + (i % 3) for i, g in enumerate(lengths)}   # 2 .. ~28, rising with the length
    units = [GeneUnit("b_0", g, n, True) for g, n in counts.items()]
    removed = {g for (_lab, g) in gene_removal(units, lengths, mode="stratified")}
    longest = sorted(lengths, key=lengths.get)[-20:]
    assert sum(g in removed for g in longest) <= 4


def test_the_other_modes_are_unchanged_by_the_stratified_arguments():
    lengths, counts = _pool()
    counts["g05"] = 60
    units = [GeneUnit("b_0", g, n, g == "g05") for g, n in counts.items()]
    base = gene_removal(units, lengths, mode="both")
    assert gene_removal(units, lengths, mode="both", strata=3, k_far=9.0, k_near=9.0) == base


def test_filter_caas_genes_in_stratified_mode(tmp_path):
    lengths, counts = _pool()
    counts.update({"g05": 60, "g06": 60})
    rows = [(g, pos) for g, n in counts.items() for pos in range(n)]
    pd.DataFrame(rows, columns=["Gene", "Position"]).to_csv(tmp_path / "discovery.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": list(lengths), "length": list(lengths.values())}).to_csv(tmp_path / "lengths.tsv", sep="\t", index=False)
    pd.DataFrame({"Gene": ["g05", "g05", "g05"], "Position": [0, 1, 2], "clustering_flag": ["Discarded"] * 3}).to_csv(
        tmp_path / "clusters.tsv", sep="\t", index=False)
    r = subprocess.run([sys.executable, str(FILTER), "-i", str(tmp_path / "discovery.tsv"), "-l", str(tmp_path / "lengths.tsv"),
                        "-c", str(tmp_path / "clusters.tsv"), "-m", "stratified", "--strata", "10", "--k-far", "3", "--k-near", "1.5",
                        "-o", str(tmp_path / "out.tsv"), "-s", str(tmp_path / "removed.tsv"), "-g", str(tmp_path / "stats.tsv")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    removed = pd.read_csv(tmp_path / "removed.tsv", sep="\t").set_index("Gene")["category"].to_dict()
    assert removed["g05"] == "Both" and removed["g06"] == "Extreme"
    stats = pd.read_csv(tmp_path / "stats.tsv", sep="\t")
    assert set(stats["category"]) <= {"Normal", "Extreme", "Dubious", "Both"}
    assert stats.loc[stats.Gene == "g05", "threshold_extreme"].iloc[0] > stats.loc[stats.Gene == "g05", "threshold_dubious"].iloc[0]
    assert not set(pd.read_csv(tmp_path / "out.tsv", sep="\t")["Gene"]) & set(removed)
