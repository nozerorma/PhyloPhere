"""Readers of perm_pos_cycle_caas.tsv.gz use the caas_score column as written and reject files without it."""
import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
ENR = ROOT / "subworkflows/ENRICHMENT/local/src"


def _write(path, rows, columns):
    with gzip.open(path, "wt", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join("" if v is None else repr(v) if isinstance(v, float) else str(v) for v in r) + "\n")


def _module(name):
    spec = importlib.util.spec_from_file_location(name, ENR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ENR))
    try:
        spec.loader.exec_module(mod)
    except ImportError as e:
        pytest.skip(f"{name} dependencies not available: {e}")
    return mod


def test_enrichment_loader_reads_caas_score_exactly_and_fills_missing_with_zero(tmp_path):
    mod = _module("posenrich_enrich")
    exact = 0.46300735781502145          # pandas' default float parser returns ...214
    _write(tmp_path / "n.tsv.gz",
           [("gA", 1, "top", "b_1", exact, 3), ("gA", 2, "top", "b_1", None, 0)],
           ["Gene", "Position", "side", "cycle", "caas_score", "n_schemes"])
    df, cycles = mod.load_caas_cycle_null(str(tmp_path / "n.tsv.gz"))
    assert df.sort_values("pos_id")["score"].tolist() == [exact, 0.0]
    assert list(cycles) == ["b_1"]


def test_enrichment_loader_rejects_a_null_without_caas_score(tmp_path):
    mod = _module("posenrich_enrich")
    _write(tmp_path / "old.tsv.gz", [("gA", 1, "top", "b_1", 0.5, 1)],
           ["Gene", "Position", "side", "cycle", "caas_sum", "n_schemes"])
    with pytest.raises(ValueError, match="caas_score"):
        mod.load_caas_cycle_null(str(tmp_path / "old.tsv.gz"))


def test_prep_script_rejects_a_null_without_caas_score(tmp_path):
    _write(tmp_path / "old.tsv.gz", [("gA", 1, "top", "b_1", 0.5, 1)],
           ["Gene", "Position", "side", "cycle", "caas_sum", "n_schemes"])
    r = subprocess.run([sys.executable, str(ENR / "posenrich_prep_caas_null.py"), "--caas-cycle-null",
                        str(tmp_path / "old.tsv.gz"), "--output", str(tmp_path / "o.pkl")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "caas_score" in r.stderr


def test_prep_script_reads_caas_score_exactly(tmp_path):
    exact = 0.46300735781502145
    _write(tmp_path / "n.tsv.gz", [("gA", 1, "top", "b_1", exact, 3)],
           ["Gene", "Position", "side", "cycle", "caas_score", "n_schemes"])
    r = subprocess.run([sys.executable, str(ENR / "posenrich_prep_caas_null.py"), "--caas-cycle-null",
                        str(tmp_path / "n.tsv.gz"), "--output", str(tmp_path / "o.pkl")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    import pickle
    d = pickle.load(open(tmp_path / "o.pkl", "rb"))
    scores = pd.concat([v["score"] for k, v in d.items() if isinstance(v, pd.DataFrame)]).tolist()
    assert exact in scores


def test_a_null_table_without_rows_is_no_null_for_the_enrichment(tmp_path):
    """N = 0 writes the header only: the loader answers 'no CAAS null' and the prep script writes an empty artifact."""
    mod = _module("posenrich_enrich")
    _write(tmp_path / "empty.tsv.gz", [], ["Gene", "Position", "side", "cycle", "caas_score", "n_schemes"])
    assert mod.load_caas_cycle_null(str(tmp_path / "empty.tsv.gz")) == (None, None)
    r = subprocess.run([sys.executable, str(ENR / "posenrich_prep_caas_null.py"), "--caas-cycle-null", str(tmp_path / "empty.tsv.gz"),
                        "--output", str(tmp_path / "o.pkl")], capture_output=True, text=True)
    assert r.returncode == 0 and "empty CAAS null" in (r.stdout + r.stderr)
    assert mod.load_prepped_caas_null(str(tmp_path / "o.pkl")) == (None, None)
