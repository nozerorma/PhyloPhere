"""core.driver.load_gene_context: cache first (canonical, then legacy layout, rst only needed), compute on a miss.

PEPC's cached ASR comes from the frozen observed inputs (golden/pepc_c4_complete/observed_inputs.tar.gz).
`run_asr_pipeline` and `codeml_slot` are replaced by spies where the question is whether PAML would run.
"""
import contextlib
import logging
import os
import shutil
import sys
import tarfile
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SRC = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1])) / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(SRC))
from src.asr.asr_single import SingleGeneASRConfig, load_precomputed_asr  # noqa: E402
from src.core import driver  # noqa: E402

GOLD = HERE / "golden/pepc_c4_complete"


@pytest.fixture(scope="module")
def pepc(tmp_path_factory):
    d = tmp_path_factory.mktemp("asr")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "align").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "align/PEPC.fasta")
    return d


def _load(pepc, cache, **kw):
    i = pepc / "observed_inputs"
    return driver.load_gene_context("PEPC", str(pepc / "align"), str(i / "pruned_tree_file.nwk"), str(i / "taxid.tsv"),
                                    "lg", str(cache), 0.1, **kw)


class Spies:
    """Stand-ins for the two PAML entry points: counts calls, returns the cached PEPC result or raises."""

    def __init__(self, monkeypatch, pepc, fail=False):
        self.runs, self.slots, self.threads = 0, 0, []
        real = self

        def run_asr_pipeline(gene, config, skip_if_exists=True, alignment_data=None, tree_data=None):
            real.runs += 1
            real.threads.append(config.threads)
            if fail:
                raise RuntimeError("codeml failed")
            cfg = SingleGeneASRConfig(alignment_path=config.alignment_path, tree_path=config.tree_path, model="lg",
                                      posterior_threshold=0.1, output_dir=pepc / "observed_inputs/asr_cache")
            return load_precomputed_asr(gene, cfg, alignment_data)

        @contextlib.contextmanager
        def codeml_slot():
            real.slots += 1
            yield

        monkeypatch.setattr(driver, "run_asr_pipeline", run_asr_pipeline)
        monkeypatch.setattr(driver, "codeml_slot", codeml_slot, raising=False)


def test_a_cached_gene_is_loaded_without_running_paml(pepc, monkeypatch):
    spies = Spies(monkeypatch, pepc)
    ctx = _load(pepc, pepc / "observed_inputs/asr_cache")
    assert ctx and ctx["node_posteriors"].posteriors_node and ctx["tree_data"].nodes
    assert (spies.runs, spies.slots) == (0, 0)


def test_the_legacy_cache_layout_is_found(pepc, monkeypatch, tmp_path):
    legacy = tmp_path / "cache/PEPC"
    shutil.copytree(pepc / "observed_inputs/asr_cache/asr_PEPC", legacy / "asr_PEPC")
    spies = Spies(monkeypatch, pepc)
    ctx = _load(pepc, tmp_path / "cache")
    assert ctx and ctx["node_posteriors"].posteriors_node
    assert spies.runs == 0


def test_a_cache_without_rst1_is_enough(pepc, monkeypatch, tmp_path):
    shutil.copytree(pepc / "observed_inputs/asr_cache/asr_PEPC", tmp_path / "cache/asr_PEPC")
    (tmp_path / "cache/asr_PEPC/rst1").unlink()
    spies = Spies(monkeypatch, pepc)
    assert _load(pepc, tmp_path / "cache")
    assert spies.runs == 0


def test_the_cached_posteriors_are_those_the_precomputed_loader_returns(pepc, monkeypatch):
    Spies(monkeypatch, pepc)
    ctx = _load(pepc, pepc / "observed_inputs/asr_cache")
    cfg = SingleGeneASRConfig(alignment_path=pepc / "align/PEPC.fasta", tree_path=pepc / "observed_inputs/pruned_tree_file.nwk",
                              model="lg", posterior_threshold=0.1, output_dir=pepc / "observed_inputs/asr_cache")
    expected = load_precomputed_asr("PEPC", cfg, ctx["alignment_data"]).posteriors_node
    assert ctx["node_posteriors"].posteriors_node == expected


def test_a_gene_missing_from_the_cache_is_computed_inside_a_codeml_slot(pepc, monkeypatch, tmp_path):
    spies = Spies(monkeypatch, pepc)
    ctx = _load(pepc, tmp_path / "empty_cache", threads=3)
    assert ctx and ctx["node_posteriors"].posteriors_node
    assert (spies.runs, spies.slots, spies.threads) == (1, 1, [3])


def test_a_gene_whose_asr_fails_is_left_out_with_a_warning(pepc, monkeypatch, tmp_path, caplog):
    spies = Spies(monkeypatch, pepc, fail=True)
    with caplog.at_level(logging.WARNING):
        assert _load(pepc, tmp_path / "empty_cache") is None
    assert spies.runs == 1 and any("PEPC" in r.getMessage() and r.levelno >= logging.WARNING for r in caplog.records)


def _job(pepc, cache, rows=()):
    i = pepc / "observed_inputs"
    # (gene, rows, alignment_dir, tree, taxid, model, cache, threshold, ensembl, trait_pairs, pss, master fields)
    return ("PEPC", list(rows), str(pepc / "align"), str(i / "pruned_tree_file.nwk"), str(i / "taxid.tsv"), "lg", str(cache), 0.1, None, {}, None, ["gene"])


def test_the_observed_step_leaves_out_a_gene_without_asr_and_does_not_fail(pepc, monkeypatch, tmp_path, caplog):
    import observed_b0_main as obs

    monkeypatch.setattr(obs, "load_gene_context", lambda *a, **k: None)
    with caplog.at_level(logging.INFO):
        assert obs._score_gene(_job(pepc, tmp_path / "cache")) == ("PEPC", None)
    assert not any(r.levelno >= logging.ERROR for r in caplog.records)


def test_the_observed_step_loads_the_asr_through_the_core_loader_with_the_cache_it_was_given(pepc, monkeypatch, tmp_path):
    import observed_b0_main as obs

    seen = []
    monkeypatch.setattr(obs, "load_gene_context", lambda *a, **k: seen.append((a, k)))
    obs._score_gene(_job(pepc, tmp_path / "cache"))
    args, _ = seen[0]
    assert args[0] == "PEPC" and args[5] == str(tmp_path / "cache") and args[6] == 0.1  # gene, cache directory, posterior threshold


def test_a_failure_in_one_gene_is_logged_and_gives_no_rows(pepc, monkeypatch, tmp_path, caplog):
    import observed_b0_main as obs

    def boom(*a, **k):
        raise RuntimeError("scoring failed")

    monkeypatch.setattr(obs, "load_gene_context", boom)
    with caplog.at_level(logging.ERROR):
        assert obs._score_gene(_job(pepc, tmp_path / "cache")) == ("PEPC", None)
    assert any("PEPC" in r.getMessage() and "scoring failed" in r.getMessage() for r in caplog.records)
