"""The observed disambiguation reproduces the frozen PEPC master CSV, with or without decoration.

golden/pepc_c4_complete/observed_inputs.tar.gz holds what the pipeline fed CT_DISAMBIGUATION_RUN for
the Tier 1 PEPC genotypic run (global_meta_caas.tsv, the 100 traitfiles and their hypothesis pairs, the
cached ASR of PEPC, tree, taxid map, gene list). This is the safety net for changes to the observed path.
"""
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GOLD = HERE / "golden/pepc_c4_complete"
MAIN = ROOT / "subworkflows/CT_DISAMBIGUATION/local/disambiguation_main.py"


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    d = tmp_path_factory.mktemp("observed")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align" / "PEPC.fasta")
    return d


def _run(inputs, out, *extra):
    i = inputs / "observed_inputs"
    cmd = [sys.executable, str(MAIN), "--alignment-dir", str(inputs / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
           "--caas-metadata", str(i / "global_meta_caas.tsv"), "--trait-file", str(i / "traitfiles"),
           "--output-dir", str(out), "--asr-mode", "compute", "--asr-model", "lg", "--posterior-threshold", "0.1",
           "--threads", "1", "--workers", "2", "--max-tasks-per-child", "50",
           "--hypotheses-pairs", str(i / "traitfiles/contrast_hypotheses_pairs.tsv"),
           "--asr-cache-dir", str(i / "asr_cache"), "--taxid-mapping", str(i / "taxid.tsv"),
           "--ensembl-genes-file", str(i / "gene_ensembl.tsv"), *extra]
    p = subprocess.run(cmd, cwd=MAIN.parent, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    return pd.read_csv(out / "caas_convergence_master.csv", keep_default_na=False)


def _golden():
    return pd.read_csv(GOLD / "caas_convergence_master.csv", keep_default_na=False)


def _assert_same_up_to_float_noise(got, gold, tol=1e-12):
    """Every cell equal, except float columns, which may differ by summation rounding (<= tol).

    The pooled scores add M hypothesis scores; the sum is correctly rounded now, so a pool of many
    hypotheses can differ from the naive sum of the frozen run in the last bits.
    """
    worst = 0.0
    for c in gold.columns:
        if gold[c].dtype.kind == "f":
            d = (got[c] - gold[c]).abs()
            assert (got[c].isna() == gold[c].isna()).all(), c
            worst = max(worst, float(d.max(skipna=True) or 0.0))
        else:
            assert got[c].equals(gold[c]), c
    print(f"max float difference to the frozen run: {worst:.2e}")
    assert worst <= tol


def test_master_csv_matches_the_frozen_pipeline_run(inputs, tmp_path):
    got = _run(inputs, tmp_path / "out")
    gold = _golden()
    assert len(gold) == 217 and list(got.columns) == list(gold.columns) and len(gold.columns) == 46
    _assert_same_up_to_float_noise(got, gold)


@pytest.mark.skipif(not os.environ.get("RUN_SLOW"), reason="decoration costs ~80 s; set RUN_SLOW=1")
def test_decoration_does_not_change_the_master_csv(inputs, tmp_path):
    _assert_same_up_to_float_noise(_run(inputs, tmp_path / "out", "--run-diagnostics", "--verbose"), _golden())
