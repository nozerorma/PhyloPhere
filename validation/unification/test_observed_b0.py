"""The observed labeling through the core: discovery rows -> full records -> master rows.

The input is the frozen PEPC discovery.tab (100 hypotheses, 5 schemes) with PEPC's cached ASR; the reference is the
frozen caas_convergence_master.csv (217 rows x 46 columns) of the pipeline run that made them. CAAS ids are content
hashes now, so `tag_support` is compared by structure (how many ids, with which counts, and that each id is the hash
of a row at that position and scheme), every other column exactly (floats within 1e-12).
"""
import csv
import gzip
import os
import random
import re
import shutil
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SRC = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1])) / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(SRC))
from src.core.driver import load_gene_context  # noqa: E402
from src.core.labelings import observed_pss, read_trait_pairs  # noqa: E402
from src.core.master import write_master_csv  # noqa: E402
from src.core.meta import caas_id  # noqa: E402
from src.core.observed import observed_entries, observed_master_rows, score_observed, unresolved_entries  # noqa: E402
from src.core.master import master_fields  # noqa: E402
from frozen_master import without_new  # noqa: E402

GOLD = HERE / "golden/pepc_c4_complete"


@pytest.fixture(scope="module")
def pepc(tmp_path_factory):
    d = tmp_path_factory.mktemp("b0")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align/PEPC.fasta")
    i = d / "observed_inputs"
    ctx = load_gene_context("PEPC", str(d / "align"), str(i / "pruned_tree_file.nwk"), str(i / "taxid.tsv"), "lg",
                            str(i / "asr_cache"), 0.1)
    rows = pd.read_csv(gzip.open(GOLD / "discovery.tab.gz", "rt"), sep="\t", dtype=str, keep_default_na=False).to_dict("records")
    return dict(dir=d, ctx=ctx, rows=rows, trait_pairs=read_trait_pairs(i / "traitfiles"),
                pss=observed_pss(i / "traitfiles/contrast_hypotheses_pairs.tsv"), fields=master_fields(4))


def _master(pepc, rows, out):
    entries = observed_entries("PEPC", rows)
    results = score_observed(pepc["ctx"], "PEPC", entries, pepc["trait_pairs"], pepc["pss"], 0.1)
    write_master_csv(observed_master_rows("PEPC", results, pepc["fields"]), out, pepc["fields"])
    return pd.read_csv(out, keep_default_na=False)


def _gold():
    return pd.read_csv(GOLD / "caas_convergence_master.csv", keep_default_na=False)


def _same_but_tag_support(got, gold, tol=1e-12):
    got = without_new(got)       # the frozen master predates the columns of frozen_master.NEW_COLUMNS
    assert list(got.columns) == list(gold.columns) and len(got) == len(gold)
    for c in gold.columns:
        if c == "tag_support":
            continue
        if gold[c].dtype.kind == "f":
            assert (got[c].isna() == gold[c].isna()).all(), c
            assert float((got[c] - gold[c]).abs().max(skipna=True) or 0.0) <= tol, c
        else:
            assert got[c].equals(gold[c]), c


def _tally(cell):
    return [(i, int(n)) for i, n in (p.rsplit(":", 1) for p in str(cell).split(",") if p)]


def test_the_b0_rows_give_the_frozen_master_apart_from_the_ids(pepc, tmp_path):
    got, gold = _master(pepc, pepc["rows"], tmp_path / "m.csv"), _gold()
    assert len(gold) == 217 and len(gold.columns) == 46
    _same_but_tag_support(got, gold)


def test_tag_support_keeps_its_shape_and_holds_the_content_ids_of_the_rows_at_that_position(pepc, tmp_path):
    got, gold = _master(pepc, pepc["rows"], tmp_path / "m.csv"), _gold()
    by_cell = {}
    for r in pepc["rows"]:
        by_cell.setdefault((int(r["position"]), r["caap_group"]), set()).add(
            caas_id("PEPC", r["position"], r["trait"], r["caap_group"], r["caas"], r["amino_encoded"], r["pattern"]))
    for g, o in zip(got.itertuples(), gold.itertuples()):
        mine, theirs = _tally(g.tag_support), _tally(o.tag_support)
        assert sorted(n for _, n in mine) == sorted(n for _, n in theirs)  # the same tally of counts
        assert {i for i, _ in mine} <= by_cell[(g.msa_pos, g.caap_group)]
        assert all(re.fullmatch(r"CAAS_[0-9A-F]{16}", i) for i, _ in mine)


def test_entries_carry_the_content_id_the_hypothesis_and_the_parsed_conserved_pair(pepc):
    first = {}
    for r in pepc["rows"]:  # one row per distinct (is_conserved_meta, conserved_pair): TRUE and FALSE, '0:', '1:2', '2:1,2', ...
        first.setdefault((r["is_conserved_meta"], r["conserved_pair"]), r)
    rows = list(first.values())
    assert {r["is_conserved_meta"] for r in rows} == {"TRUE", "FALSE"} and len(rows) >= 6
    for r, e in zip(rows, observed_entries("PEPC", rows)):
        assert e.tag == caas_id("PEPC", r["position"], r["trait"], r["caap_group"], r["caas"], r["amino_encoded"], r["pattern"])
        assert (e.position, e.position_one_based, e.trait, e.caap_group, e.caas) == (int(r["position"]), int(r["position"]) + 1,
                                                                                    r["trait"], r["caap_group"], r["caas"])
        assert e.conserved_pair == (r["conserved_pair"].split(":", 1)[-1] if r["conserved_pair"] else "")
        assert e.is_conserved_meta == (r["is_conserved_meta"] in ("TRUE", "True", "true", "1"))


def test_the_row_order_of_the_input_changes_nothing_in_the_master_ties_included(pepc, tmp_path):
    shuffled = list(pepc["rows"])
    random.Random(3).shuffle(shuffled)
    key = ["msa_pos", "caap_group", "side"]
    a = _master(pepc, pepc["rows"], tmp_path / "a.csv").sort_values(key, kind="stable").reset_index(drop=True)
    b = _master(pepc, shuffled, tmp_path / "b.csv").sort_values(key, kind="stable").reset_index(drop=True)
    assert a.equals(b)
    # every row says whether its pool had a tied derived residue (PEPC has none that decides an agreement: its tied raw
    # residues D and N share an encoding)
    assert set(a["agreement_ambiguous"].astype(str)) == {"False"}


def test_a_result_flagged_as_ambiguous_reaches_the_master_row(pepc):
    from types import SimpleNamespace
    from src.core.master import master_row
    from src.utils.gene_wrapper import convert_convergence_result_to_dict
    base = next(iter(score_observed(pepc["ctx"], "PEPC", observed_entries("PEPC", pepc["rows"][:40]), pepc["trait_pairs"], pepc["pss"], 0.1)))
    for flag in (True, False):
        flagged = SimpleNamespace(**{**vars(base), "agreement_ambiguous": flag}) if hasattr(base, "__dict__") else None
        assert flagged is not None
        row = master_row(convert_convergence_result_to_dict(flagged, multi_hypothesis=None), pepc["fields"])
        assert row["agreement_ambiguous"] == str(flag)


def test_the_discovery_order_gives_the_rows_in_the_order_of_the_frozen_master(pepc, tmp_path):
    got, gold = _master(pepc, pepc["rows"], tmp_path / "m.csv"), _gold()
    assert list(zip(got.msa_pos, got.caap_group, got.side)) == list(zip(gold.msa_pos, gold.caap_group, gold.side))


def test_entries_without_a_resolvable_hypothesis_are_rejected_in_a_multi_contrast_design(pepc):
    rows = [dict(r) for r in pepc["rows"][:3]]
    rows[1]["trait"] = "traitfile.tab"          # names no hypothesis
    rows[2]["trait"] = "traitfile_H9999.tab"    # names one the design does not have
    entries = observed_entries("PEPC", rows)
    assert [e.trait for e in unresolved_entries(entries, pepc["trait_pairs"])] == ["traitfile.tab", "traitfile_H9999.tab"]
    with pytest.raises(ValueError, match="2 entries name no hypothesis"):
        score_observed(pepc["ctx"], "PEPC", entries, pepc["trait_pairs"], pepc["pss"], 0.1)


def test_a_single_contrast_design_accepts_any_trait_name(pepc):
    rows = [dict(r, trait="traitfile.tab") for r in pepc["rows"][:3]]
    assert unresolved_entries(observed_entries("PEPC", rows), {1: [("a", "b")]}) == []
