#!/usr/bin/env python3
"""Tests for the VEP convergence_schemes gate + der_aas hardening.

Run: python -m pytest test_map_gate.py   (or: python test_map_gate.py)

map_to_cosmic.py wraps its logic in main(); map_to_primateai.py is a top-level
script. Both carry an identical `anc_der_from_descriptor` / `load_convergence_skip`
pair, so the unit tests import the cosmic copy and the primateai copy is checked
by subprocess.
"""
import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent


def _load(mod_name):
    spec = importlib.util.spec_from_file_location(mod_name, _HERE / f"{mod_name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cosmic = _load("map_to_cosmic")


# ── anc_der_from_descriptor: der_aas excludes residues present on both clades ──
def test_der_excludes_ancestral():
    # top clade carries {Y, L, V}; bottom (ancestral) carries {V}; change on top.
    anc, der = cosmic.anc_der_from_descriptor("", "Y:2,L:1,V:5", "V:11", "top")
    assert anc == {"V"}
    assert der == {"Y", "L"}                 # V dropped: it did not change


def test_der_keepall_when_subtraction_empties():
    # every derived letter is also ancestral -> keep the raw derived set (safety)
    anc, der = cosmic.anc_der_from_descriptor("", "V:3", "V:11", "top")
    assert der == {"V"}


def test_both_side_no_ancestral():
    anc, der = cosmic.anc_der_from_descriptor("", "I:2", "L:2", "both")
    assert anc == set()
    assert der == {"I", "L"}


# ── load_convergence_skip ─────────────────────────────────────────────────────
def _write(p, rows):
    p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


def test_gate_collects_empty_convergence(tmp_path):
    ps = tmp_path / "position_scores.tsv"
    _write(ps, [
        ("Gene", "Position", "convergence_schemes", "CAAS_score"),
        ("BRCA1", "96", "", "0.4"),          # genuine disagreement -> skip
        ("BRCA1", "97", "GS3,GS4", "0.5"),   # chemical convergence -> keep
        ("TP53", "10", "US", "0.6"),         # identity -> keep
        ("TP53", "11", "  ", "0.2"),         # whitespace-only -> skip
    ])
    skip = cosmic.load_convergence_skip(str(ps))
    assert skip == {("BRCA1", 96), ("TP53", 11)}


def test_gate_noop_when_absent():
    assert cosmic.load_convergence_skip(None) is None
    assert cosmic.load_convergence_skip("NO_FILE") is None
    assert cosmic.load_convergence_skip("/nonexistent/path.tsv") is None


def test_gate_noop_when_column_missing(tmp_path):
    ps = tmp_path / "ps.tsv"
    _write(ps, [("Gene", "Position", "CAAS_score"), ("BRCA1", "96", "0.4")])
    assert cosmic.load_convergence_skip(str(ps)) is None


# ── primateai.py end to end (subprocess): the 5th arg is accepted + gates ──────
def _make_inputs(tmp_path, with_skip):
    caas = tmp_path / "caas.tsv"
    _write(caas, [
        ("tag", "caas", "Gene", "Position", "side", "amino_encoded",
         "caap_group", "derived_residues", "top_residue_support", "bottom_residue_support"),
        ("t1", "Y/V", "GENE1", "5", "top", "Y>V", "US", "Y/V", "Y:2", "V:2"),
        ("t2", "F/L", "GENE2", "9", "top", "F>L", "US", "F/L", "F:2", "L:2"),
    ])
    ps = tmp_path / "ps.tsv"
    rows = [("Gene", "Position", "convergence_schemes")]
    if with_skip:
        rows += [("GENE1", "5", ""), ("GENE2", "9", "GS4")]
    else:
        rows += [("GENE1", "5", "US"), ("GENE2", "9", "GS4")]
    _write(ps, rows)
    mapdir = tmp_path / "maps"
    mapdir.mkdir(exist_ok=True)
    db = tmp_path / "pai.gz"
    with gzip.open(db, "wt") as fh:
        fh.write("chr\tpos\tref_aa\talt_aa\tscore_PAI3D\n")
    return caas, mapdir, db, ps


def _run_primateai(tmp_path, caas, mapdir, db, ps=None):
    out = tmp_path / "out.tsv"
    argv = [str(caas), str(mapdir), str(db), str(out)]
    if ps is not None:
        argv.append(str(ps))
    r = subprocess.run(
        [sys.executable, str(_HERE / "map_to_primateai.py"), *argv],
        capture_output=True, text=True, cwd=tmp_path)
    return r, (out.read_text() if out.exists() else "")


def test_primateai_accepts_optional_position_scores(tmp_path):
    caas, mapdir, db, ps = _make_inputs(tmp_path, with_skip=True)
    r, _ = _run_primateai(tmp_path, caas, mapdir, db)              # 4 args (legacy)
    assert r.returncode == 0, r.stderr
    assert "convergence gate is a no-op" in r.stderr
    r, _ = _run_primateai(tmp_path, caas, mapdir, db, ps)          # 5 args (gate)
    assert r.returncode == 0, r.stderr
    assert "1 position(s) will be skipped" in r.stderr


def test_primateai_gate_drops_target(tmp_path):
    caas, mapdir, db, ps = _make_inputs(tmp_path, with_skip=True)
    r, _ = _run_primateai(tmp_path, caas, mapdir, db, ps)
    # GENE1/5 gated out, GENE2/9 kept -> 1 unique target loaded
    assert "1 unique (Gene, Position) targets loaded" in r.stderr
    caas, mapdir, db, ps = _make_inputs(tmp_path, with_skip=False)
    r, _ = _run_primateai(tmp_path, caas, mapdir, db, ps)
    assert "2 unique (Gene, Position) targets loaded" in r.stderr


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    ok = True
    import tempfile
    for fn in fns:
        try:
            if fn.__code__.co_argcount:
                with tempfile.TemporaryDirectory() as d:
                    fn(Path(d))
            else:
                fn()
            print("PASS ", fn.__name__)
        except Exception:
            ok = False
            print("FAIL ", fn.__name__)
            traceback.print_exc()
    print("\n" + ("ALL TESTS PASSED" if ok else "SOME TESTS FAILED"))
    sys.exit(0 if ok else 1)
