"""Several workers ask for the ASR of the same gene on a cold cache: every one of them gets a complete context.

The null replays one gene in chunks on a pool of workers, and every chunk loads the gene's context. On a cold cache the
first worker runs PAML while the others read a cache entry that is still being written; a worker that fails to read it
must wait for the writer, not drop its chunk. Two tests simulate PAML (a slow writer that leaves a half-written rst for a
while), one runs the real codeml where it is installed.
"""
import contextlib
import multiprocessing as mp
import os
import shutil
import sys
import tarfile
import time
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SRC = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1])) / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(SRC))
from src.asr.asr_single import SingleGeneASRConfig, load_precomputed_asr  # noqa: E402
from src.core import driver  # noqa: E402

GOLD = HERE / "golden/pepc_c4_complete"
needs_codeml = pytest.mark.skipif(shutil.which("codeml") is None, reason="codeml not on PATH")


@pytest.fixture(scope="module")
def pepc(tmp_path_factory):
    d = tmp_path_factory.mktemp("asrc")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align/PEPC.fasta")
    return d


def _load(pepc, cache):
    i = pepc / "observed_inputs"
    return driver.load_gene_context("PEPC", str(pepc / "align"), str(i / "pruned_tree_file.nwk"), str(i / "taxid.tsv"),
                                    "lg", str(cache), 0.1)


def _slow_writer(pepc, delay=1.5):
    """A stand-in for run_asr_pipeline: leaves a half-written rst for `delay` seconds, then the complete entry."""
    source = pepc / "observed_inputs/asr_cache/asr_PEPC"

    def run_asr_pipeline(gene, config, skip_if_exists=True, alignment_data=None, tree_data=None):
        out = Path(config.output_dir) / f"asr_{gene}"
        out.mkdir(parents=True, exist_ok=True)
        data = (source / "rst").read_bytes()
        (out / "rst").write_bytes(data[: len(data) // 2])
        time.sleep(delay)
        shutil.copytree(source, out, dirs_exist_ok=True)
        cfg = SingleGeneASRConfig(alignment_path=config.alignment_path, tree_path=config.tree_path, model="lg",
                                  posterior_threshold=0.1, output_dir=Path(config.output_dir))
        return load_precomputed_asr(gene, cfg, alignment_data)

    return run_asr_pipeline


def _worker(args):
    pepc, cache, barrier_wait = args
    time.sleep(barrier_wait)
    ctx = _load(Path(pepc), cache)
    return ctx is not None and bool(ctx["node_posteriors"].posteriors_node)


def _four_workers(pepc, cache):
    ctx = mp.get_context("fork")           # the children inherit the patched module
    with ctx.Pool(4) as pool:
        return pool.map(_worker, [(str(pepc), str(cache), 0.05 * k) for k in range(4)], chunksize=1)


def test_workers_that_arrive_during_the_computation_wait_for_it(pepc, tmp_path, monkeypatch):
    monkeypatch.setattr(driver, "run_asr_pipeline", _slow_writer(pepc))
    assert _four_workers(pepc, tmp_path / "cold_cache") == [True] * 4


def test_the_asr_is_computed_once_for_the_four_workers(pepc, tmp_path, monkeypatch):
    runs = tmp_path / "runs.txt"
    inner = _slow_writer(pepc)

    def counted(*a, **kw):
        with open(runs, "a") as fh:
            fh.write("x\n")
        return inner(*a, **kw)

    monkeypatch.setattr(driver, "run_asr_pipeline", counted)
    assert _four_workers(pepc, tmp_path / "cold_cache") == [True] * 4
    assert runs.read_text().count("x") == 1


def test_a_truncated_entry_left_by_a_crashed_run_is_recomputed(pepc, tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    shutil.copytree(pepc / "observed_inputs/asr_cache/asr_PEPC", cache / "asr_PEPC")
    rst = cache / "asr_PEPC/rst"
    rst.write_bytes(rst.read_bytes()[:2000])
    seen = []
    inner = _slow_writer(pepc, delay=0)

    def spy(gene, config, skip_if_exists=True, alignment_data=None, tree_data=None):
        seen.append(skip_if_exists)
        return inner(gene, config, skip_if_exists, alignment_data, tree_data)

    monkeypatch.setattr(driver, "run_asr_pipeline", spy)
    ctx = _load(pepc, cache)
    assert ctx and ctx["node_posteriors"].posteriors_node
    assert seen == [False]                 # a half-written entry must not be taken for a finished one


@needs_codeml
def test_four_workers_on_a_cold_cache_with_the_real_codeml(pepc, tmp_path):
    assert _four_workers(pepc, tmp_path / "cold_cache") == [True] * 4
