"""core.master: the master CSV is written from the workers' rows, in the order the database export used."""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.core.master import master_row, serialize_value, write_master_csv  # noqa: E402

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


def test_process_all_genes_needs_max_pairs():
    from src.utils.gene_wrapper import process_all_genes
    with pytest.raises(ValueError, match="max_pairs"):
        process_all_genes(genes=["G"], alignment_dir="a", tree_file="t", caas_metadata_path="m", trait_file_path="f",
                          taxid_mapping_path=None, asr_mode="compute", asr_model="lg", asr_cache_dir=None,
                          posterior_threshold=0.1, threads_per_gene=1, workers=1, run_diagnostics=False,
                          output_dir=Path("."), max_pairs=None)
