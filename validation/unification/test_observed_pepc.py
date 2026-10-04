"""The observed labeling end to end on PEPC: discovery rows -> master shard -> ct_disambiguation/caas_convergence_master.csv.

golden/pepc_c4_complete/observed_inputs.tar.gz holds what the pipeline fed the observed scoring for the Tier 1 PEPC
genotypic run (the 100 traitfiles and their hypothesis pairs, the cached ASR of PEPC, tree, taxid map, gene list) and
discovery.tab.gz the 8 205 discovery rows; caas_convergence_master.csv is the frozen master of that run. The ids of
`tag_support` are content hashes now (test_observed_b0 checks their shape); every other column equals the frozen one,
floats within 1e-12 (a pool of many hypotheses is summed with math.fsum, which can differ from the frozen naive sum
in the last bits).
"""
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
GOLD = HERE / "golden/pepc_c4_complete"
LOCAL = ROOT / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(HERE))
from frozen_master import without_new  # noqa: E402
from test_observed_b0 import _gold, _same_but_tag_support  # noqa: E402


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    d = tmp_path_factory.mktemp("observed")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align" / "PEPC.fasta")
    with gzip.open(GOLD / "discovery.tab.gz", "rt") as src, open(d / "discovery.tab", "w") as dst:
        shutil.copyfileobj(src, dst)
    return d


def _observed(inputs, out, workers):
    """CAAS_OBSERVED's two steps by hand: score the discovery rows, then write the master."""
    i = inputs / "observed_inputs"
    shards = out / "shards"
    steps = [
        [sys.executable, str(LOCAL / "observed_b0_main.py"), "--alignment-dir", str(inputs / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
         "--discovery", str(inputs / "discovery.tab"), "--design", str(i / "traitfiles"), "--output-dir", str(shards), "--asr-model", "lg",
         "--posterior-threshold", "0.1", "--workers", str(workers), "--asr-cache-dir", str(i / "asr_cache"), "--taxid-mapping", str(i / "taxid.tsv"),
         "--ensembl-genes-file", str(i / "gene_ensembl.tsv"), "--fop-pairs", str(i / "traitfiles/contrast_hypotheses_pairs.tsv")],
        [sys.executable, str(LOCAL / "contract_main.py"), "--b0-dirs", str(shards), "--design", str(i / "traitfiles"),
         "--discovery-file", str(inputs / "discovery.tab"), "--output-dir", str(out)],
    ]
    for cmd in steps:
        p = subprocess.run(cmd, capture_output=True, text=True)
        assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    return out / "ct_disambiguation/caas_convergence_master.csv"


def test_master_csv_matches_the_frozen_pipeline_run(inputs, tmp_path):
    got = pd.read_csv(_observed(inputs, tmp_path, 2), keep_default_na=False)
    gold = _gold()
    assert len(gold) == 217 and list(without_new(got).columns) == list(gold.columns) and len(gold.columns) == 46
    _same_but_tag_support(got, gold)


def test_the_master_does_not_depend_on_the_number_of_workers(inputs, tmp_path):
    one = _observed(inputs, tmp_path / "one", 1).read_bytes()
    two = _observed(inputs, tmp_path / "two", 2).read_bytes()
    assert one == two and len(one) > 100_000
