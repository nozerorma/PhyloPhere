"""core.master: the master CSV is written from the workers' rows, in the order the database export used."""
import csv
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.labelings import design_max_pairs  # noqa: E402
from src.core.master import master_fields, master_row, serialize_value, write_master_csv  # noqa: E402

FIELDS = ["gene", "msa_pos", "side", "asr_path_score", "participating_hypotheses"]


def _rows(path):
    return list(csv.DictReader(open(path, newline="")))


def test_serialization_matches_the_database_export():
    assert serialize_value(None) == "" and serialize_value(True) == "True" and serialize_value(0.5) == "0.5"
    assert serialize_value(["H1", "H2"]) == "H1,H2" and serialize_value(("a", 3)) == "a,3"
    assert master_row({"gene": "G", "side": None, "extra": 1}, FIELDS) == {
        "gene": "G", "msa_pos": "", "side": "", "asr_path_score": "", "participating_hypotheses": ""}


def test_rows_are_ordered_by_gene_then_position_and_keep_production_order(tmp_path):
    def r(gene, pos, side):
        return (gene, pos, master_row({"gene": gene, "msa_pos": pos, "side": side}, FIELDS))
    rows = [r("B", 10, "top"), r("A", 9, "bottom"), r("A", 100, "top"), r("A", 9, "top"), r("B", None, "none")]
    n = write_master_csv(rows, tmp_path / "m.csv", FIELDS)
    got = _rows(tmp_path / "m.csv")
    assert n == 5
    assert [(x["gene"], x["msa_pos"], x["side"]) for x in got] == [
        ("A", "9", "bottom"), ("A", "9", "top"),  # same position: as the worker produced them
        ("A", "100", "top"),                       # numeric, not lexical: 9 before 100
        ("B", "", "none"), ("B", "10", "top")]     # a missing position sorts first, as NULL does in SQL


def test_empty_input_writes_only_the_header(tmp_path):
    assert write_master_csv([], tmp_path / "m.csv", FIELDS) == 0
    assert open(tmp_path / "m.csv").read().strip() == ",".join(FIELDS)


def test_the_master_columns_follow_the_number_of_pairs_of_the_design():
    one, four = master_fields(1), master_fields(4)
    assert four[:len(one)] == one and len(four) - len(one) == 3 * 8 and len(four) == 46
    assert [f for f in four if f.startswith("domain_4_")] == [
        f"domain_4_{k}" for k in ("posterior", "score", "anc_aa", "top_aa", "bot_aa", "anc_aa_support", "top_aa_support", "bot_aa_support")]
    golden = (Path(__file__).resolve().parent / "golden/pepc_c4_complete/caas_convergence_master.csv").read_text().splitlines()[0]
    assert ",".join(four) == golden  # the frozen master has the columns of a four-pair design


def test_the_pair_count_of_a_design_is_its_largest_pair_id(tmp_path):
    (tmp_path / "dir").mkdir()
    (tmp_path / "dir/traitfile_H1.tab").write_text("a\t1\t1\nb\t0\t1\nc\t1\t3\nd\t0\t3\n")
    (tmp_path / "dir/traitfile_H2.tab").write_text("a\t1\t2\nb\t0\t2\n")
    (tmp_path / "one.tab").write_text("a\t1\t1\nb\t0\t1\ne\t1\t2\nf\t0\t2\n")
    (tmp_path / "bad.tab").write_text("a\t1\nb\t0\n")
    assert design_max_pairs(tmp_path / "dir") == 3 and design_max_pairs(tmp_path / "one.tab") == 2
    assert design_max_pairs(tmp_path / "bad.tab") == 1 and design_max_pairs(tmp_path / "missing") == 1  # nothing readable: one pair


def test_a_gz_path_is_written_compressed_with_the_same_content(tmp_path):
    rows = [("A", 9, master_row({"gene": "A", "msa_pos": 9, "side": "top"}, FIELDS)),
            ("A", 2, master_row({"gene": "A", "msa_pos": 2, "side": "bottom"}, FIELDS))]
    write_master_csv(rows, tmp_path / "m.csv", FIELDS)
    write_master_csv(rows, tmp_path / "m.csv.gz", FIELDS)
    assert gzip.open(tmp_path / "m.csv.gz", "rt").read() == open(tmp_path / "m.csv").read()
