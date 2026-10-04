"""compare_contract.py: the observed contract files of two runs, by set of rows, with ids recomputed and ties counted.

Synthetic runs for the table rules, the frozen PEPC master for the master rules, and the whole PEPC chain (frozen
files and the former report's meta tables as run A; the b_0 path with shuffled discovery rows as run B).
"""
import csv
import gzip
import io
import os
import random
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
LOCAL = ROOT / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LOCAL))
import compare_contract as cc  # noqa: E402
from src.core import contract  # noqa: E402
from test_contract import _rmd_meta, needs_r  # noqa: E402

GOLD = HERE / "golden/pepc_c4_complete"
DISC_TEXT = gzip.open(GOLD / "discovery.tab.gz", "rt").read()
MASTER_TEXT = (GOLD / "caas_convergence_master.csv").read_text()
H = "gene\tmode\tcaap_group\ttrait\tposition\tcaas\tamino_encoded\tpattern\tffgn\tfbgn\tgfg\tgbg\tmfg\tmbg\tffg\tfbg\tms"
ROWS = ["A1\tCAAP\tUS\ttraitfile_H1.tab\t3\tA/B\tA/B\t1\t2\t2\t0\t0\t0\t0\ts1,s2\ts3,s4\ts9,s8",
        "A1\tCAAP\tGS1\ttraitfile_H1.tab\t3\tA/B\tp/q\t1\t2\t2\t0\t0\t0\t0\ts1,s2\ts3,s4\ts9,s8",
        "B2\tCAAP\tUS\ttraitfile_H2.tab\t7\tC/D\tC/D\t2\t2\t2\t0\t0\t0\t0\ts1,s2\ts3,s4\t"]


def _write(run, rel, text):
    p = run / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def _tables(run, rows=ROWS, background="A1\t3,4,5\nB2\t7\n", genes="A1\nB2\n"):
    _write(run, cc.DISCOVERY, "\n".join([H] + rows) + "\n")
    _write(run, cc.BACKGROUND, background)
    _write(run, cc.BACKGROUND_GENES, genes)
    return run


def _status(results, rel):
    return results[rel]["pass"]


# ── discovery, background, genes ─────────────────────────────────────────────

def test_identical_runs_pass_and_files_absent_from_both_are_skipped(tmp_path):
    r = cc.compare_runs(_tables(tmp_path / "a"), _tables(tmp_path / "b"))
    assert [_status(r, k) for k in (cc.DISCOVERY, cc.BACKGROUND, cc.BACKGROUND_GENES)] == [True] * 3
    assert _status(r, cc.MASTER) is None


def test_the_row_order_and_the_order_of_the_missing_species_do_not_matter(tmp_path):
    swapped = [ROWS[1], ROWS[0].replace("s9,s8", "s8,s9"), ROWS[2]]  # same genes in order, rows swapped inside a gene, ms reordered
    r = cc.compare_runs(_tables(tmp_path / "a"), _tables(tmp_path / "b", rows=swapped))
    assert _status(r, cc.DISCOVERY) is True


@pytest.mark.parametrize("rows", [ROWS[:2], ROWS[:2] + [ROWS[2].replace("\t7\t", "\t8\t")], ROWS + [ROWS[0]]])
def test_a_missing_a_changed_or_an_extra_row_fails(tmp_path, rows):
    r = cc.compare_runs(_tables(tmp_path / "a"), _tables(tmp_path / "b", rows=rows))
    assert _status(r, cc.DISCOVERY) is False


def test_b_must_be_ordered_by_gene_but_a_is_not_required_to_be(tmp_path):
    reordered = [ROWS[2], ROWS[0], ROWS[1]]
    assert _status(cc.compare_runs(_tables(tmp_path / "a"), _tables(tmp_path / "b", rows=reordered)), cc.DISCOVERY) is False
    assert _status(cc.compare_runs(_tables(tmp_path / "c", rows=reordered), _tables(tmp_path / "d")), cc.DISCOVERY) is True


def test_a_different_header_or_a_file_missing_from_one_run_fails(tmp_path):
    a = _tables(tmp_path / "a")
    b = _tables(tmp_path / "b")
    _write(b, cc.DISCOVERY, "\n".join([H.replace("\tms", "\tms2")] + ROWS) + "\n")
    assert _status(cc.compare_runs(a, b), cc.DISCOVERY) is False
    (b / cc.BACKGROUND).unlink()
    r = cc.compare_runs(a, b)
    assert r[cc.BACKGROUND] == {"pass": False, "missing_from": "B"}


def test_background_positions_and_gene_lists_are_compared(tmp_path):
    a = _tables(tmp_path / "a")
    assert _status(cc.compare_runs(a, _tables(tmp_path / "b", background="A1\t3,4\nB2\t7\n")), cc.BACKGROUND) is False
    assert _status(cc.compare_runs(a, _tables(tmp_path / "c", background="B2\t7\nA1\t3,4,5\n")), cc.BACKGROUND) is False  # B unordered
    assert _status(cc.compare_runs(a, _tables(tmp_path / "d", genes="A1\n")), cc.BACKGROUND_GENES) is False


# ── master ───────────────────────────────────────────────────────────────────

def _master_run(run, text):
    return _write(run, cc.MASTER, text) and run


def _frame(text=MASTER_TEXT):
    return pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)


def _text(df):
    return df.to_csv(index=False)


def _tie_cell(df):
    """(row, column, the other tied residue) of a modal-residue cell whose support has two residues at the maximum."""
    for c in df.columns:
        if re.fullmatch(r"domain_\d+_(anc|top|bot)_aa", c):
            for i, cell in df[c + "_support"].items():
                t = cc._tally(cell)
                top = [k for k, n in t.items() if n == max(t.values())] if t else []
                if len(top) >= 2 and df.at[i, c] in top:
                    return i, c, next(k for k in top if k != df.at[i, c]), t
    raise AssertionError("no tie in the frozen master")


def _master_result(tmp_path, mutate):
    df = _frame()
    mutate(df)
    r = cc.compare_runs(_master_run(tmp_path / "a", MASTER_TEXT), _master_run(tmp_path / "b", _text(df)))
    return r[cc.MASTER]


def test_the_frozen_master_equals_itself_in_any_row_order(tmp_path):
    df = _frame().sample(frac=1, random_state=1)
    r = cc.compare_runs(_master_run(tmp_path / "a", MASTER_TEXT), _master_run(tmp_path / "b", _text(df)))
    assert r[cc.MASTER]["pass"] is True and r[cc.MASTER]["n_a"] == 217


def test_floats_pass_within_the_tolerance_and_fail_above_it(tmp_path):
    def shift(delta):
        def f(df):
            df["asr_path_score"] = [repr(float(x) + delta) if x else x for x in df["asr_path_score"]]
        return f
    assert _master_result(tmp_path / "ok", shift(1e-14))["pass"] is True
    r = _master_result(tmp_path / "bad", shift(1e-9))
    assert r["pass"] is False and r["worst_column"] == "asr_path_score" and r["max_abs_delta"] > 1e-12


def test_a_modal_residue_that_differs_where_two_residues_tie_is_tolerated_and_counted(tmp_path):
    i, c, other, _ = _tie_cell(_frame())
    r = _master_result(tmp_path, lambda df: df.at.__setitem__((i, c), other))
    assert r["pass"] is True and r["modal_residue_tie_cells_tolerated"] == 1


def test_a_modal_residue_that_differs_without_a_tie_fails(tmp_path):
    i, c, _, tally = _tie_cell(_frame())
    r = _master_result(tmp_path, lambda df: df.at.__setitem__((i, c), "@"))
    assert r["pass"] is False and r["modal_residue_cells_differing_without_a_tie"] == 1


def test_any_other_text_column_a_missing_row_and_a_changed_tally_shape_fail(tmp_path):
    assert _master_result(tmp_path / "1", lambda df: df.at.__setitem__((0, "caas"), "Z/Z"))["pass"] is False
    assert _master_result(tmp_path / "2", lambda df: df.drop(index=3, inplace=True))["pass"] is False
    def drop_one_id(df):
        first = df.at[0, "tag_support"].split(",")
        df.at[0, "tag_support"] = ",".join(first[1:])
    r = _master_result(tmp_path / "3", drop_one_id)
    assert r["pass"] is False and r["tag_support_shape_differs"] == 1


def _renamed_ids(df):
    """The same ids under other names: one new id per old id, so the tallies keep their shape."""
    import hashlib
    return df.assign(tag_support=df["tag_support"].map(lambda t: re.sub(r"CAAS_[0-9A-Z]+", lambda m: "CAAS_" + hashlib.sha256(m.group(0).encode()).hexdigest()[:16].upper(), t)))


def test_a_different_id_of_the_same_shape_passes_without_meta_and_fails_when_it_is_not_in_the_meta_of_b(tmp_path):
    assert _master_result(tmp_path / "plain", lambda df: df.update(_renamed_ids(df)))["pass"] is True  # nothing to check the ids against
    a = _master_run(tmp_path / "a", MASTER_TEXT)
    b = _master_run(tmp_path / "b", _text(_renamed_ids(_frame())))
    (tmp_path / "d.tab").write_text(DISC_TEXT)
    contract.write_meta(tmp_path / "d.tab", b / cc.META_DIR)
    r = cc.compare_runs(a, b)[cc.MASTER]
    assert r["pass"] is False and r["tag_support_ids_outside_meta_b"] > 0


# ── meta tables ──────────────────────────────────────────────────────────────

def _meta_runs(tmp_path):
    (tmp_path / "d.tab").write_text(DISC_TEXT)
    b = tmp_path / "b"
    contract.write_meta(tmp_path / "d.tab", b / cc.META_DIR)
    a = tmp_path / "a"
    shutil.copytree(b / cc.META_DIR, a / cc.META_DIR)
    for f in (a / cc.META_DIR).glob("*.tsv"):  # the former report's ids were seeded draws, different strings
        df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        df["tag"] = [f"CAAS_AAAAA{i:04d}A" for i in range(len(df))]
        df.to_csv(f, sep="\t", index=False)
    return a, b


def test_meta_tables_are_compared_without_the_tag_and_every_new_tag_is_recomputed(tmp_path):
    a, b = _meta_runs(tmp_path)
    r = cc.compare_runs(a, b)
    names = [k for k in r if k.startswith(cc.META_DIR)]
    assert len(names) == 6 and all(r[k]["pass"] is True for k in names) and r[f"{cc.META_DIR}/global_meta_caas.tsv"]["ids_checked"] == 8205


def test_a_tag_that_is_not_the_content_id_a_missing_file_and_a_changed_row_fail(tmp_path):
    a, b = _meta_runs(tmp_path)
    g = b / cc.META_DIR / "global_meta_caas.tsv"
    df = pd.read_csv(g, sep="\t", dtype=str, keep_default_na=False)
    df.at[0, "tag"] = "CAAS_" + "0" * 16
    df.to_csv(g, sep="\t", index=False)
    r = cc.compare_runs(a, b)[f"{cc.META_DIR}/global_meta_caas.tsv"]
    assert r["pass"] is False and r["ids_not_the_content_id"] == 1
    (b / cc.META_DIR / "US_meta_caas.tsv").unlink()
    assert cc.compare_runs(a, b)[f"{cc.META_DIR}/US_meta_caas.tsv"] == {"pass": False, "missing_from": "B"}
    df2 = pd.read_csv(a / cc.META_DIR / "GS1_meta_caas.tsv", sep="\t", dtype=str, keep_default_na=False)
    df2.at[0, "Position"] = "99999"
    df2.to_csv(a / cc.META_DIR / "GS1_meta_caas.tsv", sep="\t", index=False)
    assert cc.compare_runs(a, b)[f"{cc.META_DIR}/GS1_meta_caas.tsv"]["pass"] is False


# ── the command line ─────────────────────────────────────────────────────────

def test_the_command_line_prints_a_line_per_file_writes_the_report_and_sets_the_exit_status(tmp_path):
    a, b = _tables(tmp_path / "a"), _tables(tmp_path / "b")
    cmd = [sys.executable, str(HERE / "compare_contract.py"), "--a", str(a), "--b", str(b), "--report", str(tmp_path / "r.json")]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0 and "[PASS] caastools/discovery.tab" in p.stdout and "3 passed" in p.stdout
    assert (tmp_path / "r.json").exists()
    _write(b, cc.BACKGROUND, "A1\t3\nB2\t7\n")
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 1 and "[FAIL] caastools/background.output" in p.stdout


# ── the whole PEPC chain ─────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def pepc_runs(tmp_path_factory):
    d = tmp_path_factory.mktemp("chain")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align/PEPC.fasta")
    i = d / "observed_inputs"
    # run A: the frozen files and the meta tables of the former report
    a = d / "a"
    _write(a, cc.DISCOVERY, DISC_TEXT)
    _write(a, cc.BACKGROUND, gzip.open(GOLD / "background.output.gz", "rt").read())
    genes = sorted({l.split("\t")[0] for l in (a / cc.BACKGROUND).read_text().splitlines() if l.strip()})
    _write(a, cc.BACKGROUND_GENES, "".join(g + "\n" for g in genes))
    _write(a, cc.MASTER, MASTER_TEXT)
    for name, text in _rmd_meta(d / "rmd", DISC_TEXT).items():
        _write(a, f"{cc.META_DIR}/{name}", text)
    # run B: shuffled discovery rows (the order of the entries differs) through the b_0 path
    lines = DISC_TEXT.splitlines()
    body = lines[1:]
    random.Random(3).shuffle(body)
    shuffled = d / "shuffled.tab"
    shuffled.write_text("\n".join([lines[0]] + body) + "\n")
    shards, b = d / "shards", d / "b"
    steps = [[sys.executable, str(LOCAL / "observed_b0_main.py"), "--alignment-dir", str(d / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
              "--discovery", str(shuffled), "--design", str(i / "traitfiles"), "--output-dir", str(shards), "--asr-model", "lg",
              "--posterior-threshold", "0.1", "--workers", "2", "--asr-cache-dir", str(i / "asr_cache"), "--taxid-mapping", str(i / "taxid.tsv"),
              "--ensembl-genes-file", str(i / "gene_ensembl.tsv"), "--fop-pairs", str(i / "traitfiles/contrast_hypotheses_pairs.tsv")],
             [sys.executable, str(LOCAL / "contract_main.py"), "--b0-dirs", str(shards), "--design", str(i / "traitfiles"),
              "--discovery-file", str(shuffled), "--output-dir", str(b)]]
    for cmd in steps:
        p = subprocess.run(cmd, capture_output=True, text=True)
        assert p.returncode == 0, p.stdout[-1000:] + p.stderr[-1000:]
    shutil.move(str(b / "meta_caas"), str(d / "meta_tmp"))  # the published layout: meta_caas/meta_caas/
    shutil.move(str(d / "meta_tmp"), str(b / cc.META_DIR))
    _write(b, cc.DISCOVERY, shuffled.read_text())
    shutil.copy(a / cc.BACKGROUND, b / cc.BACKGROUND)
    shutil.copy(a / cc.BACKGROUND_GENES, b / cc.BACKGROUND_GENES)
    return a, b


@needs_r
def test_the_pepc_chain_passes_with_shuffled_discovery_rows_and_reports_the_ties(pepc_runs):
    a, b = pepc_runs
    r = cc.compare_runs(a, b)
    assert {k: v["pass"] for k, v in r.items()} == {k: True for k in r}, {k: v for k, v in r.items() if v["pass"] is not True}
    assert len([k for k in r if k.startswith(cc.META_DIR)]) == 6
    m = r[cc.MASTER]
    assert m["n_a"] == m["n_b"] == 217 and m["max_abs_delta"] <= 1e-12 and m["modal_residue_tie_cells_tolerated"] >= 1
    assert m["tag_support_ids_outside_meta_b"] == 0


@needs_r
def test_the_pepc_chain_fails_when_one_cell_of_b_is_wrong(pepc_runs, tmp_path):
    a, b = pepc_runs
    shutil.copytree(b, tmp_path / "b2")
    f = tmp_path / "b2" / cc.MASTER
    df = pd.read_csv(f, dtype=str, keep_default_na=False)
    df.at[10, "asr_path_score"] = repr(float(df.at[10, "asr_path_score"]) + 1e-6)
    df.to_csv(f, index=False)
    r = cc.compare_runs(a, tmp_path / "b2")
    assert r[cc.MASTER]["pass"] is False and all(v["pass"] is True for k, v in r.items() if k != cc.MASTER)
