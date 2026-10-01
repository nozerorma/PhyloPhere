"""compare_null.py: identifiers exact, values within the tolerance, the RDS by all.equal, missing tables detected."""
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_null as cn  # noqa: E402

TOY = HERE / "cancer_b0_toy/neoplasia_prevalence_toy_complete"


def _run(tmp, scores, quantiles="cycle\tscheme\tq50\nb_1\tUS\t0.5\nb_1\tGS1\t0.25\n"):
    d = tmp / "caas_permulation"
    d.mkdir(parents=True)
    (d / "gene_cycle_scores.tsv").write_text("Gene\tcycle\tgene_caas_score\tn\n" + "".join(f"{g}\t{c}\t{v!r}\t{n}\n" for g, c, v, n in scores))
    (d / "perm_pos_quantiles.tsv").write_text(quantiles)
    return tmp


SCORES = [("A", "b_1", 0.1, 3), ("A", "b_2", 0.30000000000000004, 3), ("B", "b_1", 1.0, 2), ("B", "b_2", 0.0, 2)]


def _results(a, b, **kw):
    return cn.compare_runs(a, b, **kw)


def test_identical_runs_pass_and_absent_tables_are_skipped(tmp_path):
    r = _results(_run(tmp_path / "a", SCORES), _run(tmp_path / "b", SCORES))
    assert r["caas_permulation/gene_cycle_scores.tsv"]["pass"] is True
    assert r["caas_permulation/perm_pos_sample.tsv"]["pass"] is None


def test_a_value_difference_above_the_tolerance_fails_and_below_it_passes(tmp_path):
    a = _run(tmp_path / "a", SCORES)
    below = _run(tmp_path / "below", [("A", "b_1", 0.1 + 1e-14, 3)] + SCORES[1:])
    above = _run(tmp_path / "above", [("A", "b_1", 0.1 + 1e-6, 3)] + SCORES[1:])
    ok = _results(a, below)["caas_permulation/gene_cycle_scores.tsv"]
    bad = _results(a, above)["caas_permulation/gene_cycle_scores.tsv"]
    assert ok["pass"] is True and ok["n_value_cells_not_bitwise_equal"] == 1 and ok["max_abs_delta"] < 1e-12
    assert bad["pass"] is False and bad["worst_column"] == "gene_caas_score" and bad["max_abs_delta"] > 9e-7


def test_row_order_does_not_matter_but_a_missing_or_extra_row_does(tmp_path):
    a = _run(tmp_path / "a", SCORES)
    assert _results(a, _run(tmp_path / "shuffled", SCORES[::-1]))["caas_permulation/gene_cycle_scores.tsv"]["pass"] is True
    missing = _results(a, _run(tmp_path / "missing", SCORES[:-1]))["caas_permulation/gene_cycle_scores.tsv"]
    assert missing["pass"] is False and missing["only_a"] == 1
    extra = _results(a, _run(tmp_path / "extra", SCORES + [("C", "b_1", 0.5, 1)]))["caas_permulation/gene_cycle_scores.tsv"]
    assert extra["pass"] is False and extra["only_b"] == 1


def test_an_identifier_difference_fails_even_when_values_agree(tmp_path):
    other = [("A", "b_1", 0.1, 3), ("A", "b_3", 0.30000000000000004, 3)] + SCORES[2:]
    r = _results(_run(tmp_path / "a", SCORES), _run(tmp_path / "b", other))["caas_permulation/gene_cycle_scores.tsv"]
    assert r["pass"] is False and r["only_a"] == 1 and r["only_b"] == 1


def test_an_na_in_one_run_and_a_number_in_the_other_fails(tmp_path):
    a = _run(tmp_path / "a", SCORES)
    b = _run(tmp_path / "b", SCORES)
    p = b / "caas_permulation/gene_cycle_scores.tsv"
    p.write_text(p.read_text().replace("b_1\t1.0\t2", "b_1\tNA\t2"))
    r = _results(a, b)["caas_permulation/gene_cycle_scores.tsv"]
    assert r["pass"] is False and r["na_mismatch"] == 1


def test_a_table_missing_from_one_run_fails(tmp_path):
    a = _run(tmp_path / "a", SCORES)
    b = _run(tmp_path / "b", SCORES)
    (b / "caas_permulation/perm_pos_quantiles.tsv").unlink()
    assert _results(a, b)["caas_permulation/perm_pos_quantiles.tsv"] == {"pass": False, "missing_from": "B"}


def test_the_b0_slice_and_extra_tables_are_compared(tmp_path):
    a, b = _run(tmp_path / "a", SCORES), _run(tmp_path / "b", SCORES)
    for run, v in ((a, 0.5), (b, 0.5 + 1e-6)):
        (run / "caas_permulation/b0").mkdir()
        (run / "caas_permulation/b0/gene_cycle_scores.tsv").write_text(f"Gene\tcycle\tgene_caas_score\nA\tb_0\t{v!r}\n")
        (run / "scoring").mkdir()
        (run / "scoring/gene_scores.tsv").write_text("gene\tp\nA\t0.5\n")
    r = _results(a, b, extra=["scoring/gene_scores.tsv"])
    assert r["caas_permulation/b0/gene_cycle_scores.tsv"]["pass"] is False and r["scoring/gene_scores.tsv"]["pass"] is True


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not on PATH")
def test_the_rds_is_compared_with_all_equal(tmp_path):
    a, b = _run(tmp_path / "a", SCORES), _run(tmp_path / "b", SCORES)
    for run, v in ((a, 1.0), (b, 1.0 + 1e-14)):
        subprocess.run(["Rscript", "-e", f'saveRDS(list(m = matrix(c({v!r}, 2, 3, 4), 2), gene_stat = "size_adj_max"), "{run}/caas_permulation/caas_perms.rds")'], check=True)
    assert _results(a, b)["caas_permulation/caas_perms.rds"]["pass"] is True
    subprocess.run(["Rscript", "-e", f'saveRDS(list(m = matrix(c(1.1, 2, 3, 4), 2), gene_stat = "size_adj_max"), "{b}/caas_permulation/caas_perms.rds")'], check=True)
    assert _results(a, b)["caas_permulation/caas_perms.rds"]["pass"] is False


@pytest.mark.skipif(not TOY.exists(), reason="cancer_b0_toy fixture not present")
def test_a_stored_run_equals_itself_and_a_perturbed_copy_differs(tmp_path):
    run = tmp_path / "copy"
    shutil.copytree(TOY / "caas_permulation", run / "caas_permulation", ignore=shutil.ignore_patterns("perm_disc"), symlinks=False)
    same = cn.compare_runs(TOY, run, rscript="Rscript-not-installed-here")
    assert all(r["pass"] in (True, None) for r in same.values()) and same["caas_permulation/gene_cycle_scores.tsv"]["n_a"] > 100
    p = run / "caas_permulation/gene_cycle_scores.tsv"
    lines = p.read_text().splitlines()
    cols = lines[0].split("\t")
    i = cols.index("global_asr")
    row = lines[1].split("\t")
    row[i] = repr(float(row[i]) + 1e-6)
    p.write_text("\n".join([lines[0], "\t".join(row)] + lines[2:]) + "\n")
    assert cn.compare_runs(TOY, run, rscript="x")["caas_permulation/gene_cycle_scores.tsv"]["pass"] is False
