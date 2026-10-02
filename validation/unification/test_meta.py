"""core.meta: a CAAS id is a function of the content of one discovery.tab row, and of nothing else."""
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SRC = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1])) / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(SRC))
from src.core.meta import ID_HEX_CHARS, assign_ids, caas_id  # noqa: E402

TOY = HERE / "cancer_b0_toy/neoplasia_prevalence_toy_complete/caastools/discovery.tab"
GOLD = HERE / "golden/pepc_c4_complete/discovery.tab.gz"
ROW = dict(gene="C3orf38", position=80, hypothesis="H10", caap_group="US", caas="IIVV/IIII", amino_encoded="IIVV/IIII", pattern="3")


def test_the_id_has_a_fixed_shape():
    assert re.fullmatch(rf"CAAS_[0-9A-F]{{{ID_HEX_CHARS}}}", caas_id(**ROW)) and ID_HEX_CHARS == 16


def test_the_same_content_gives_the_same_id_in_any_process():
    here = caas_id(**ROW)
    code = ("import sys; sys.path.insert(0, %r); from src.core.meta import caas_id; "
            "print(caas_id(%s))" % (str(SRC), ", ".join(f"{k}={v!r}" for k, v in ROW.items())))
    for seed in ("1", "2", "random"):  # str.__hash__ is randomized per process; the id must not follow it
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED=seed))
        assert out.stdout.strip() == here, out.stderr


@pytest.mark.parametrize("field, other", [("gene", "C3orf39"), ("position", 81), ("hypothesis", "H11"), ("caap_group", "GS1"),
                                          ("caas", "IIVV/IIIV"), ("amino_encoded", "hhtt/hhhh"), ("pattern", "2")])
def test_every_field_of_the_row_is_part_of_the_id(field, other):
    assert caas_id(**{**ROW, field: other}) != caas_id(**ROW)


def test_the_hypothesis_is_normalized_and_the_position_is_an_integer():
    ids = {caas_id(**{**ROW, "hypothesis": h}) for h in ("H10", "traitfile_H10.tab", "b_0~H10")}
    assert len(ids) == 1
    assert caas_id(**{**ROW, "position": "80"}) == caas_id(**ROW)


def test_fields_cannot_run_into_each_other():
    # adjacent text fields: concatenated they would read the same
    assert caas_id("G", 1, "H1", "US", "AB", "C", "1") != caas_id("G", 1, "H1", "US", "A", "BC", "1")
    assert caas_id("G", 1, "H1", "US", "AB/CD", "x", "1") != caas_id("G", 1, "H1", "US", "AB", "/CDx", "1")


def test_an_id_does_not_depend_on_the_other_rows():
    rows = [tuple({**ROW, "position": p}.values()) for p in range(50)]
    full = assign_ids(rows)
    assert assign_ids(rows[::-1])[::-1] == full  # order
    assert assign_ids(rows[10:20]) == full[10:20]  # which other rows are present
    assert full[7] == caas_id(**{**ROW, "position": 7})


def test_two_different_rows_with_the_same_id_are_an_error():
    rows = [tuple({**ROW, "position": p}.values()) for p in range(400)]
    with pytest.raises(ValueError, match="collision"):
        assign_ids(rows, hex_chars=2)  # 256 possible ids for 400 distinct rows


def test_identical_rows_share_their_id_without_error():
    row = tuple(ROW.values())
    assert len(set(assign_ids([row, row, row]))) == 1


def _ids_of(df):
    df = df.copy()
    df["hyp"] = df["trait"]
    return assign_ids(zip(df["gene"], df["position"], df["hyp"], df["caap_group"], df["caas"], df["amino_encoded"], df["pattern"]))


@pytest.mark.skipif(not TOY.exists(), reason="cancer_b0_toy fixture not present")
def test_the_ids_of_a_stored_cancer_toy_run_are_unique():
    df = pd.read_csv(TOY, sep="\t", dtype=str, keep_default_na=False)
    ids = _ids_of(df)
    assert len(df) > 3000 and len(set(ids)) == len(df)


def test_the_ids_of_the_pepc_golden_discovery_are_unique():
    df = pd.read_csv(gzip.open(GOLD, "rt"), sep="\t", dtype=str, keep_default_na=False)
    ids = _ids_of(df)
    assert len(df) > 8000 and len(set(ids)) == len(df)
