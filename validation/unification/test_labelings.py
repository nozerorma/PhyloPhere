"""core.labelings: one reader of the labelings, the observed design and the PSS weights."""
import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core import labelings as L  # noqa: E402

GOLD = Path(__file__).resolve().parent / "golden/pepc_c4_complete"
B0_SCRIPT = ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"


def test_ids_and_bases():
    assert L.hyp_id("b_12~H3") == "H3" and L.hyp_id("traitfile_H10.tab") == "H10" and L.hyp_id("b_7") == "H1"
    assert L.base_cycle("b_12~H3") == "b_12" and L.base_cycle("b_12") == "b_12"
    lab = L.Labeling("b_4~H2", ("f1", "f2"), ("b1", "b2"))
    assert (lab.base, lab.hyp) == ("b_4", "H2")
    assert lab.trait_pairs() == {1: [("f1", "b1"), ("f2", "b2")]} == L.trait_pairs_from(lab.fg, lab.bg)
    assert L.Labeling("b_4", ("f",), ("b",)).hyp == "H1"


def test_read_labelings_prefers_fop_file(tmp_path):
    (tmp_path / "resample_001.tab").write_text("b_1\ta,b\tc,d\n")
    assert list(L.read_labelings(tmp_path)) == ["b_1"]
    (tmp_path / "fop_labelings.tab").write_text("b_1~H1\ta,b\tc,d\nb_1~H2\tb,a\td,c\nbad\tonly\n")
    got = L.read_labelings(tmp_path)
    assert list(got) == ["b_1~H1", "b_1~H2"] and got["b_1~H2"].fg == ("b", "a")
    assert L.read_labelings(tmp_path / "resample_001.tab")["b_1"].bg == ("c", "d")  # a single file works too


def _traitfile(path, pairs):
    path.write_text("".join(f"{f}\t1\t{k}\n{b}\t0\t{k}\n" for k, (f, b) in pairs))


def test_read_design_orders_pairs_numerically_and_hypotheses_by_file_name(tmp_path):
    _traitfile(tmp_path / "traitfile_H10.tab", [(2, ("f2", "b2")), (1, ("f1", "b1"))])  # pair ids out of order
    _traitfile(tmp_path / "traitfile_H2.tab", [(1, ("g1", "c1"))])
    (tmp_path / "traitfile_fop.tab").write_text("x\t1\t1\ny\t0\t1\n")  # not a hypothesis file
    got = L.read_design(tmp_path)
    # hypotheses in file-name order (H10 before H2), as the observed disambiguation has always visited them:
    # the order enters floating-point sums downstream
    assert list(got) == ["b_0~H10", "b_0~H2"]
    assert got["b_0~H10"].fg == ("f1", "f2") and got["b_0~H10"].bg == ("b1", "b2")  # pairs by pair id
    single = L.read_design(tmp_path / "traitfile_H2.tab")
    assert list(single) == ["b_0"] and single["b_0"].fg == ("g1",)


def test_read_trait_pairs_groups_the_pairs_by_hypothesis_in_file_name_order(tmp_path):
    _traitfile(tmp_path / "traitfile_H1.tab", [(1, ("f1", "b1")), (2, ("f2", "b2"))])
    _traitfile(tmp_path / "traitfile_H3.tab", [(1, ("g1", "c1"))])
    assert L.read_trait_pairs(tmp_path) == {1: [("f1", "b1"), ("f2", "b2")], 3: [("g1", "c1")]}


def test_read_pss_both_formats_and_sentinels(tmp_path):
    fop = tmp_path / "fop_pairs.tsv"
    fop.write_text("cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\nb_1\tH1\t1\ta\tb\t2.5\nb_1\tH2\t1\ta\tb\tNA\nb_2\tH1\t2\ta\tb\t1\n")
    got = L.read_pss(fop)
    assert got == {"b_1": {("H1", 1): 2.5}, "b_2": {("H1", 2): 1.0}}  # the NA weight is dropped
    obs = tmp_path / "contrast_hypotheses_pairs.tsv"
    obs.write_text("hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\ntraitfile_H3\t2\ta\tb\t4\n")
    assert L.observed_pss(obs) == {("H3", 2): 4.0}  # the id is normalized, the cycle defaults to b_0
    assert L.read_pss("NO_FOP_PAIRS") == {} and L.read_pss(tmp_path / "missing.tsv") == {}
    assert L.observed_pss("NO_HYP_PAIRS") is None
    (tmp_path / "other.tsv").write_text("a\tb\n1\t2\n")
    assert L.read_pss(tmp_path / "other.tsv") == {}


def test_golden_pepc_pss_and_b0_script_agree_with_the_reader(tmp_path):
    pss = L.observed_pss(GOLD / "contrast_hypotheses_pairs.tsv")
    assert len(pss) == 400 and {h for h, _d in pss} == {f"H{i}" for i in range(1, 101)}
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (cfg / "contrast_hypotheses_pairs.tsv").write_text((GOLD / "contrast_hypotheses_pairs.tsv").read_text())
    by_hyp = {}
    for r in csv.DictReader(open(cfg / "contrast_hypotheses_pairs.tsv"), delimiter="\t"):
        by_hyp.setdefault(r["hypothesis_id"], []).append(r)
    for h, rows in by_hyp.items():
        _traitfile(cfg / f"traitfile_{h}.tab", [(int(r["pair"]), (r["species1"], r["species2"])) for r in rows])
    design = L.read_design(cfg)
    assert len(design) == 100
    subprocess.run([sys.executable, str(B0_SCRIPT), "--config", str(cfg), "--fop", "--labelings-out", str(tmp_path / "b0.tab"),
                    "--pairs-out", str(tmp_path / "pairs.tsv")], check=True, capture_output=True)
    assert L.read_labelings(tmp_path / "b0.tab") == design  # the script and the reader see the same design
    # its fop_pairs rows carry the same weights as the observed file, under the cycle b_0
    (tmp_path / "hdr.tsv").write_text("cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n" + (tmp_path / "pairs.tsv").read_text())
    assert L.read_pss(tmp_path / "hdr.tsv") == {"b_0": pss}
    subprocess.run([sys.executable, str(B0_SCRIPT), "--config", str(cfg), "--labelings-out", str(tmp_path / "plain.tab")],
                   check=True, capture_output=True)
    assert list(L.read_labelings(tmp_path / "plain.tab")) == ["b_0"]
    assert L.read_labelings(tmp_path / "plain.tab")["b_0"].fg == design["b_0~H1"].fg
