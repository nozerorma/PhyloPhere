"""trains_coordinates.py: trains in trimmed versus untrimmed alignment coordinates."""
import random
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
import trains_coordinates as tc  # noqa: E402
from src.core.postproc import ctrain  # noqa: E402


def _map(ori_cols):
    """prot_ali_col (1-based) -> untrimmed column for the trimmed columns 1..n."""
    return {i + 1: c for i, c in enumerate(ori_cols)}


def _unit(positions, ori, **kw):
    return tc.unit_flags(positions, _map(ori), **kw)


def test_a_triplet_that_is_adjacent_only_because_columns_were_removed_is_not_a_train_untrimmed():
    # trimmed columns 10, 11, 12 (0-based 9..11) sit at untrimmed columns 10, 14, 20
    ori = list(range(1, 10)) + [10, 14, 20]
    comps, ft, fu = _unit([9, 10, 11], ori)
    assert ft == {9, 10, 11} and fu == set()
    assert comps == [[9, 10, 11]]


def test_one_removed_column_inside_a_triplet_keeps_it_two_do_not():
    # untrimmed 1, 2, 4: 3 of 4 = 0.75 >= 0.7 stays a train; untrimmed 1, 3, 5: 3 of 5 = 0.6 does not
    assert _unit([0, 1, 2], [1, 2, 4])[2] == {0, 1, 2}
    assert _unit([0, 1, 2], [1, 3, 5])[2] == set()


def test_without_removed_columns_both_coordinates_agree():
    pos = [3, 4, 5, 9, 10, 11, 12, 30]
    ori = list(range(1, 41))
    comps, ft, fu = _unit(pos, ori)
    assert ft == fu == {3, 4, 5, 9, 10, 11, 12}


def test_components_hold_exactly_the_positions_ctrain_flags():
    rng = random.Random(7)
    for _ in range(200):
        pos = rng.sample(range(60), rng.randint(0, 25))
        maxcaas = rng.choice([0.5, 0.7, 0.9])
        minlen = rng.choice([2, 3, 4])
        comps = tc.train_components(pos, maxcaas, minlen)
        assert sorted(p for c in comps for p in c) == ctrain(pos, maxcaas, minlen)
        assert all(len(c) >= 2 for c in comps)
        assert len({p for c in comps for p in c}) == sum(len(c) for c in comps)  # components are disjoint


def test_components_split_at_a_gap_and_merge_when_intervals_overlap():
    assert tc.train_components([1, 2, 3, 20, 21, 22]) == [[1, 2, 3], [20, 21, 22]]
    assert tc.train_components([1, 2, 3, 4, 5]) == [[1, 2, 3, 4, 5]]


def test_untrimmed_flags_are_a_subset_of_trimmed_flags_with_the_default_parameters():
    rng = random.Random(11)
    for _ in range(300):
        n = 80
        removed = set(rng.sample(range(1, 121), 40))
        ori = [c for c in range(1, 121) if c not in removed][:n]
        pos = rng.sample(range(len(ori)), rng.randint(0, 30))
        _, ft, fu = _unit(pos, ori)
        assert fu <= ft


def test_with_a_lower_maxcaas_untrimmed_coordinates_can_flag_what_trimmed_ones_do_not():
    # two CAAS with one removed column between: trimmed span 2 < minlen, untrimmed 2 of 3 = 0.67 >= 0.6;
    # a distant third position keeps the unit above ctrain's n >= minlen entry condition
    ori = [1, 3] + list(range(10, 59))
    _, ft, fu = _unit([0, 1, 50], ori, maxcaas=0.6, minlen=3)
    assert ft == set() and fu == {0, 1}


def test_measure_counts_components_by_size_and_outcome_and_skips_genes_without_a_map():
    rows = [("A", "US", "H1", p) for p in (9, 10, 11)]                      # compressed triplet: lost
    rows += [("A", "GS1", "H1", p) for p in (0, 1, 2, 3)]                   # quartet without removed columns: whole
    rows += [("B", "US", "H1", p) for p in (0, 1, 2)]                       # no map
    df = pd.DataFrame(rows, columns=["gene", "caap_group", "trait", "position"])
    ori_a = list(range(1, 10)) + [10, 14, 20]
    comps, tot, skipped = tc.measure(df, {"A": _map(ori_a)})
    got = {(r.caap_group, r.bin): r.status for r in comps.itertuples()}
    assert got == {("US", "3"): "lost", ("GS1", "4"): "whole"}
    assert tot["flagged_trimmed"] == 7 and tot["flagged_untrimmed"] == 4 and tot["only_untrimmed"] == 0
    assert tot["records_removed_trimmed"] == 7 and tot["records_removed_untrimmed"] == 4
    assert tot["units_with_train_trimmed"] == 2 and tot["units_with_train_untrimmed"] == 1
    assert dict(skipped) == {"B": 1}


def test_command_line_finds_a_map_named_for_another_reference_species(tmp_path, capsys, monkeypatch):
    (tmp_path / "G.Lemur_catta.map.tsv").write_text(
        "ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n" + "".join(f"{c}\tselected\t{c}\t{c}\n" for c in range(1, 8)))
    (tmp_path / "d.tsv").write_text("gene\tcaap_group\ttrait\tposition\n" + "".join(f"G\tUS\tH1\t{p}\n" for p in (1, 2, 3)))
    monkeypatch.setattr(sys, "argv", ["x", "--discovery", str(tmp_path / "d.tsv"), "--map-dir", str(tmp_path)])
    tc.main()
    out = capsys.readouterr().out
    assert "genes with a map: 1" in out and "flagged positions: trimmed 3  untrimmed 3" in out
