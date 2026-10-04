"""`ct perm-replay --export_b0_discovery`: the discovery.tab rows of the b_0 labelings, from the kernel.

The kernel decides which (position, hypothesis, scheme) are CAAS; the rows carry the 17 or 19 columns of
the scalar discovery (reference_discovery.py), named by the real traitfiles of the design and in a fixed
order (position, then trait by file name, then scheme). The scalar's own row order inside a position follows the
glob order of the trait directory, which is the filesystem's, so both sides are sorted by that key before they are
compared. `ms` (missing species) is compared as a set: the scalar builds it from a Python set, so its order changes
with the hash seed of the process, and the export writes it in pair order instead.
"""
import csv
import gzip
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
sys.path.insert(0, str(HERE))
from test_kernel_fuzz import _case  # noqa: E402

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
CT = ROOT / "subworkflows/CT/local/ct"
REF = Path(__file__).resolve().parent / "reference_discovery.py"  # the scalar discovery
GOLD = HERE / "golden/pepc_c4_complete"
B0_SCRIPT = ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"


def _run(cmd, cwd, env=None, ok=True):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True, env=dict(os.environ, **(env or {})))
    if ok:
        assert p.returncode == 0, p.stdout[-800:] + p.stderr[-800:]
    return p


def _table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


SCHEME_ORDER = ["US", "GS1", "GS2", "GS3", "GS4"]


def _canonical(df):
    """Rows in the export's order: position (numeric), trait by file name, scheme in the order of SCHEMES."""
    key = pd.DataFrame({"p": df["position"].astype(int), "t": df["trait"], "s": df["caap_group"].map(SCHEME_ORDER.index)})
    return df.loc[key.sort_values(["p", "t", "s"], kind="stable").index].reset_index(drop=True)


def _assert_same_rows(got, want, ctx=""):
    want = _canonical(want)
    assert list(got.columns) == list(want.columns), ctx
    assert len(got) == len(want), f"{len(got)} rows against {len(want)} {ctx}"
    for c in want.columns:
        if c == "ms":
            same = [set(a.split(",")) == set(b.split(",")) for a, b in zip(got[c], want[c])]
            assert all(same), f"ms differs at row {same.index(False)} {ctx}"
        else:
            assert got[c].tolist() == want[c].tolist(), f"column {c} differs {ctx}"


@pytest.fixture(scope="module")
def pepc(tmp_path_factory):
    d = tmp_path_factory.mktemp("b0d")
    cfg = d / "cfg"
    cfg.mkdir()
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    for f in (d / "observed_inputs/traitfiles").glob("traitfile_H*.tab"):
        shutil.copy(f, cfg / f.name)
    shutil.copy(GOLD / "contrast_hypotheses_pairs.tsv", cfg / "contrast_hypotheses_pairs.tsv")
    shutil.copy(GOLD / "PEPC.fasta", d / "PEPC.fasta")
    _run(["python3", B0_SCRIPT, "--config", cfg, "--fop", "--labelings-out", "b0.tab", "--pairs-out", "/dev/null"], d)
    return d


PEPC_ARGS = ["--patterns", "1,2,3", "--miss_pair", "--caap_mode", "--max_conserved", "2", "--max_fg_gaps", "0", "--max_bg_gaps", "0",
             "--max_gaps", "0", "--max_fg_miss", "0", "--max_bg_miss", "0", "--max_miss", "0"]


def test_the_pepc_b0_rows_equal_the_frozen_discovery_on_all_19_columns(pepc):
    _run([CT, "perm-replay", "-a", "PEPC.fasta", "-t", "cfg", "-s", "b0.tab", "--fmt", "fasta", *PEPC_ARGS,
          "--export_b0_discovery", "b0.discovery"], pepc)
    got = _table(pepc / "b0.discovery")
    want = pd.read_csv(gzip.open(GOLD / "discovery.tab.gz", "rt"), sep="\t", dtype=str, keep_default_na=False)
    assert len(want) > 8000 and list(want.columns)[:3] == ["gene", "mode", "caap_group"]
    _assert_same_rows(got, want)
    assert got.equals(_canonical(got))  # the export is written in the fixed order


def test_only_the_b0_labelings_are_exported_when_permuted_labelings_are_in_the_file(pepc):
    b0 = (pepc / "b0.tab").read_text()
    (pepc / "mixed.tab").write_text(b0 + b0.replace("b_0~", "b_1~") + b0.replace("b_0~", "b_2~"))
    _run([CT, "perm-replay", "-a", "PEPC.fasta", "-t", "cfg", "-s", "b0.tab", "--fmt", "fasta", *PEPC_ARGS,
          "--export_b0_discovery", "alone.discovery"], pepc)
    _run([CT, "perm-replay", "-a", "PEPC.fasta", "-t", "cfg", "-s", "mixed.tab", "--fmt", "fasta", *PEPC_ARGS,
          "--export_b0_discovery", "m.discovery", "--export_perm_discovery", "m.disc"], pepc)
    assert {c.split("~")[0] for c in _table(pepc / "m.disc")["cycle"]} == {"b_0", "b_1", "b_2"}  # the null labelings do hit
    mixed, alone = _table(pepc / "m.discovery"), _table(pepc / "alone.discovery")
    assert len(alone) > 8000 and len(mixed) == len(alone)
    _assert_same_rows(mixed, alone)


@pytest.mark.parametrize("seed", range(60))
def test_the_b0_rows_equal_the_scalar_discovery_on_random_input(tmp_path, seed):
    rng = random.Random(seed)
    thr, flags = _case(rng, tmp_path)
    _run([sys.executable, str(REF), "-a", "G1.fasta", "-t", "cfg", "-o", "s.out", "--fmt", "fasta", *flags, *thr], tmp_path)
    _run([CT, "perm-replay", "-a", "G1.fasta", "-t", "cfg", "-s", "b0.tab", "--fmt", "fasta", *flags, *thr,
          "--export_b0_discovery", "k.discovery"], tmp_path)
    scalar = tmp_path / "s.out"
    if not scalar.is_file():
        assert not (tmp_path / "k.discovery").exists()  # like the scalar, a file only when there is a row
        return
    _assert_same_rows(_table(tmp_path / "k.discovery"), _table(scalar), f"seed={seed} flags={flags} thr={thr}")


def test_the_missing_species_column_does_not_depend_on_the_hash_seed(tmp_path):
    """Six species of one hypothesis are absent from the alignment, so `ms` has six entries on every row."""
    species = [f"s{i}" for i in range(12)]
    present = species[:6]
    rng = random.Random(1)
    seqs = {s: [] for s in present}
    for _ in range(30):
        r0, r1 = rng.sample("ADKEG", 2)
        for i, s in enumerate(present):
            seqs[s].append(r0 if i % 2 == 0 else r1)
    (tmp_path / "G1.fasta").write_text("".join(f">{s}\n{''.join(q)}\n" for s, q in seqs.items()))
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    # 6 pairs: the first three are in the alignment, the last three are absent
    fg, bg = species[0:12:2], species[1:12:2]
    (cfg / "traitfile_H1.tab").write_text("".join(f"{f}\t1\t{k + 1}\n{b}\t0\t{k + 1}\n" for k, (f, b) in enumerate(zip(fg, bg))))
    (tmp_path / "b0.tab").write_text(f"b_0~H1\t{','.join(fg)}\t{','.join(bg)}\n")
    args = ["--patterns", "1,2,3", "--caap_mode", "--max_conserved", "0", "--max_fg_gaps", "9", "--max_bg_gaps", "9", "--max_gaps", "9",
            "--max_fg_miss", "9", "--max_bg_miss", "9", "--max_miss", "9"]
    outs = []
    for hs in ("1", "2", "3", "4", "5"):
        _run([CT, "perm-replay", "-a", "G1.fasta", "-t", "cfg", "-s", "b0.tab", "--fmt", "fasta", *args,
              "--export_b0_discovery", f"k{hs}.discovery"], tmp_path, env={"PYTHONHASHSEED": hs})
        outs.append((tmp_path / f"k{hs}.discovery").read_bytes())
    rows = [line.split("\t") for line in outs[0].decode().splitlines()[1:]]
    assert rows and all(len(r[16].split(",")) == 6 for r in rows)
    assert len(set(outs)) == 1  # five hash seeds, one file
    assert rows[0][16].split(",")[:3] == ["s6", "s8", "s10"]  # foreground first, in pair order


def test_the_export_needs_a_resample_file(pepc):
    p = _run([CT, "perm-replay", "-a", "PEPC.fasta", "-t", "cfg", "-s", "cfg", "--fmt", "fasta", *PEPC_ARGS,
              "--export_b0_discovery", "d.discovery"], pepc, ok=False)
    assert p.returncode != 0
