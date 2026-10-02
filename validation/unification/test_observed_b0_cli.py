"""observed_b0_main.py: the b_0 discovery rows of a perm-replay batch become one master shard per gene.

Input is the frozen PEPC discovery.tab (as `PEPC.b0.discovery.tsv`) with PEPC's cached ASR; the reference is the
frozen caas_convergence_master.csv. The PSS weights come from a fop_pairs.tsv that also holds another cycle with
different weights, so only the b_0 rows may be used.
"""
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
MAIN = ROOT / "subworkflows/CT_DISAMBIGUATION/local/observed_b0_main.py"
GOLD = HERE / "golden/pepc_c4_complete"
sys.path.insert(0, str(HERE))
from test_observed_b0 import _gold, _same_but_tag_support  # noqa: E402


@pytest.fixture(scope="module")
def inp(tmp_path_factory):
    d = tmp_path_factory.mktemp("obs_cli")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align/PEPC.fasta")
    (d / "b0").mkdir()
    with gzip.open(GOLD / "discovery.tab.gz", "rt") as src, open(d / "b0/PEPC.b0.discovery.tsv", "w") as dst:
        shutil.copyfileobj(src, dst)
    pairs = list(csv.DictReader(open(d / "observed_inputs/traitfiles/contrast_hypotheses_pairs.tsv"), delimiter="\t"))
    with open(d / "fop_pairs.tsv", "w") as fh:
        fh.write("cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n")
        for cycle, weight in (("b_0", lambda w: w), ("b_1", lambda w: 1.0 / w)):  # b_1 weights are not proportional to b_0's
            for r in pairs:
                fh.write("\t".join([cycle, r["hypothesis_id"], r["pair"], r["species1"], r["species2"], str(weight(float(r["pss_score"])))]) + "\n")
    return d


def _run(inp, out, *extra, b0="b0"):
    i = inp / "observed_inputs"
    cmd = [sys.executable, str(MAIN), "--alignment-dir", str(inp / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
           "--b0-dir", str(inp / b0), "--design", str(i / "traitfiles"), "--output-dir", str(out), "--asr-model", "lg",
           "--posterior-threshold", "0.1", "--workers", "2", "--asr-cache-dir", str(i / "asr_cache"),
           "--taxid-mapping", str(i / "taxid.tsv"), "--ensembl-genes-file", str(i / "gene_ensembl.tsv"), *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_the_shard_of_a_gene_is_its_part_of_the_frozen_master(inp, tmp_path):
    p = _run(inp, tmp_path / "out", "--fop-pairs", str(inp / "fop_pairs.tsv"))
    assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    assert [f.name for f in (tmp_path / "out").iterdir()] == ["PEPC.master.csv.gz"]
    got = pd.read_csv(tmp_path / "out/PEPC.master.csv.gz", keep_default_na=False)
    _same_but_tag_support(got, _gold())


def test_a_gene_without_an_alignment_is_left_out_and_reported(inp, tmp_path):
    rows = (inp / "b0/PEPC.b0.discovery.tsv").read_text().splitlines()
    (tmp_path / "b0").mkdir()
    (tmp_path / "b0/GHOST.b0.discovery.tsv").write_text("\n".join([rows[0]] + [r.replace("PEPC\t", "GHOST\t", 1) for r in rows[1:20]]) + "\n")
    shutil.copy(inp / "b0/PEPC.b0.discovery.tsv", tmp_path / "b0")
    (tmp_path / "ens.tsv").write_text((inp / "observed_inputs/gene_ensembl.tsv").read_text() + "GHOST\tchr1\t1\t2\t+\t970\tP04711\n")
    i = inp / "observed_inputs"
    p = subprocess.run([sys.executable, str(MAIN), "--alignment-dir", str(inp / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
                        "--b0-dir", str(tmp_path / "b0"), "--design", str(i / "traitfiles"), "--output-dir", str(tmp_path / "out"),
                        "--posterior-threshold", "0.1", "--workers", "2", "--asr-cache-dir", str(i / "asr_cache"),
                        "--taxid-mapping", str(i / "taxid.tsv"), "--ensembl-genes-file", str(tmp_path / "ens.tsv"),
                        "--fop-pairs", str(inp / "fop_pairs.tsv")], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    assert [f.name for f in (tmp_path / "out").iterdir()] == ["PEPC.master.csv.gz"]
    assert "1 left out" in p.stderr + p.stdout and "GHOST" in p.stderr + p.stdout


def test_a_directory_without_b0_hits_writes_nothing(inp, tmp_path):
    (tmp_path / "empty").mkdir()
    p = _run(inp, tmp_path / "out", b0=str(tmp_path / "empty"))
    assert p.returncode == 0, p.stderr[-800:]
    assert list((tmp_path / "out").iterdir()) == []


def test_the_ids_of_the_master_are_those_of_the_meta_tables_at_the_same_position_and_scheme(inp, tmp_path):
    sys.path.insert(0, str(ROOT / "subworkflows/CT_DISAMBIGUATION/local"))
    from src.core import contract
    p = _run(inp, tmp_path / "out", "--fop-pairs", str(inp / "fop_pairs.tsv"))
    assert p.returncode == 0, p.stderr[-800:]
    contract.write_meta(inp / "b0/PEPC.b0.discovery.tsv", tmp_path / "meta")
    ids = {}
    for line in (tmp_path / "meta/global_meta_caas.tsv").read_text().splitlines()[1:]:
        c = line.split("\t")
        ids.setdefault((int(c[3]), c[5]), set()).add(c[0])
    master = pd.read_csv(tmp_path / "out/PEPC.master.csv.gz", keep_default_na=False)
    checked = 0
    for r in master.itertuples():
        for item in str(r.tag_support).split(","):
            assert item.rsplit(":", 1)[0] in ids[(r.msa_pos, r.caap_group)]
            checked += 1
    assert checked > 200
