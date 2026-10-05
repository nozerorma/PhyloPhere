"""compare_pepc_runs.py: equal runs compare equal; a different null, a missing position or a lost cycle does not."""
import importlib.util
import os
from pathlib import Path

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
spec = importlib.util.spec_from_file_location("cpr", ROOT / "validation/tier1/scripts/compare_pepc_runs.py")
cpr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cpr)

COLS = "Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\tp.adj_sam\n"
ROWS = [("PEPC", 779, 1.0, "top", 0.001, 0.0265, 0.0003), ("PEPC", 664, 1.0, "top", 0.001, 0.0265, 0.0003), ("PEPC", 700, 0.4, "top", 0.1, 0.3, 0.3),
        ("PEPC", 700, 0.5, "bottom", 0.08, 0.3, 0.3)]


def _results(tmp_path, tag, rows=ROWS, cycles=1000):
    root = tmp_path / tag
    for trait in cpr.TRAITS:
        (root / f"{trait}_complete/scoring").mkdir(parents=True)
        (root / f"{trait}_complete/caas_permulation").mkdir(parents=True)
        (root / f"{trait}_complete/scoring/position_scores.tsv").write_text(COLS + "".join("\t".join(map(str, r)) + "\n" for r in rows))
        (root / f"{trait}_complete/caas_permulation/gene_cycle_scores.tsv").write_text("Gene\tcycle\n" + "".join(f"PEPC\tb_{i}\n" for i in range(cycles)))
    return root


def test_identical_runs_are_equal_and_a_position_is_the_best_of_its_sides(tmp_path):
    a, b = _results(tmp_path, "a"), _results(tmp_path, "b")
    r = cpr.compare(a, b, "c4")
    assert r["equal"] and r["n_shared"] == 3 and r["a"]["rows"] == 4 and r["a"]["positions"] == 3
    d, u = cpr.positions(a, "c4")
    assert u.loc[700, "score"] == 0.5 and u.loc[700, "p"] == 0.08           # best side, smallest p.emp


def test_a_different_empirical_p_is_not_equal(tmp_path):
    a = _results(tmp_path, "a")
    b = _results(tmp_path, "b", rows=[ROWS[0], ROWS[1], ("PEPC", 700, 0.4, "top", 0.2, 0.3, 0.3), ("PEPC", 700, 0.5, "bottom", 0.15, 0.3, 0.3)])
    r = cpr.compare(a, b, "c4")
    assert not r["equal"] and r["max_dp"] > 0.01 and r["p_emp_equal"] == 2


def test_a_missing_position_is_not_equal(tmp_path):
    a, b = _results(tmp_path, "a"), _results(tmp_path, "b", rows=ROWS[:2])
    r = cpr.compare(a, b, "c4")
    assert not r["equal"] and r["only_a"] == [700] and r["only_b"] == []


def test_a_null_that_lost_cycles_is_not_equal_even_when_the_positions_agree(tmp_path):
    a, b = _results(tmp_path, "a"), _results(tmp_path, "b", cycles=920)
    r = cpr.compare(a, b, "c4")
    assert not r["equal"] and (r["cycles_a"], r["cycles_b"]) == (1000, 920) and r["max_dp"] == 0
