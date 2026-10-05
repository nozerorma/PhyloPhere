"""The null replays a gene in chunks of cycles; a gene whose chunks fail only in part must stop the run.

A gene left out entirely is left out of the null as it is left out of the observed. A gene with some chunks lost keeps a null
made of part of its cycles while N still counts all of them, which biases every p.emp low without a sound. The ASR replay
workers are replaced by synthetic ones (module-level, so the pool pickles them by reference); the driver loop is the real code.
"""
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
from src.utils import gene_wrapper as gw  # noqa: E402

CYCLES = ["b_1", "b_2", "b_3"]
FAILING = {}      # gene -> cycle tags whose chunk fails


def replay(args):
    # the pool starts its workers with forkserver, so what a worker must know travels in its arguments: the labelings
    gene, tags, labelings = args[0], args[7], args[8]
    if any(f"FAIL:{gene}" in labelings[t][0] for t in tags):
        return gene, None
    return gene, [(c, []) for c in tags]


def finalize(gene, pooled, n_cycles_total, postproc_filter, minlen, maxcaas, columns):
    rows = [{"Gene": gene, "cycle": c, "Position": 10 + i, "caap_group": "GS1", "asr_path_score": 0.5, "n_detected": 1, "clust": 0, "side": "top"}
            for i, (c, _) in enumerate(pooled)]
    return gene, rows


@pytest.fixture
def stubbed(monkeypatch):
    FAILING.clear()
    monkeypatch.setattr(gw, "build_cycle_inputs", lambda disc, resample, cycles: (
        list(CYCLES), {c: ([f"FAIL:{g}" for g, tags in FAILING.items() if c in tags], []) for c in CYCLES}))
    monkeypatch.setattr(gw, "_perms_worker_replay_wrapper", replay)
    monkeypatch.setattr(gw, "_perms_worker_finalize", finalize)
    yield
    FAILING.clear()


def _run(out):
    return gw.process_all_genes_perms(
        genes=["GENEA", "GENEB"], alignment_dir="x", tree_file="x", perm_discovery_file="x", resample_dir="x", taxid_mapping_path=None,
        asr_model="lg", asr_cache_dir="x", posterior_threshold=0.0, workers=1, output_dir=out, detail_only=True,
        gene_sizes={"GENEA": 100, "GENEB": 100}, chunk_threshold=10, chunk_target_size=1)       # three chunks of one cycle per gene


def _genes(out):
    return {p.name.split(".")[0] for p in (out / "perm_pos_detail").glob("*.tsv.gz")}


def test_a_gene_with_some_chunks_lost_stops_the_run(stubbed, tmp_path):
    FAILING["GENEA"] = {"b_2"}
    with pytest.raises(RuntimeError, match=r"GENEA: 1 of 3 chunks"):
        _run(tmp_path / "o")


def test_a_gene_with_every_chunk_lost_is_left_out_and_the_rest_is_written(stubbed, tmp_path):
    FAILING["GENEA"] = set(CYCLES)
    out = _run(tmp_path / "o")
    assert _genes(out) == {"GENEB"}
    rows = gzip.open(next((out / "perm_pos_detail").glob("GENEB*.tsv.gz")), "rt").read().splitlines()
    assert len(rows) == 1 + len(CYCLES)


def test_no_failure_writes_every_gene(stubbed, tmp_path):
    assert _genes(_run(tmp_path / "o")) == {"GENEA", "GENEB"}


def test_the_consistency_rule_alone():
    gw._require_consistent_chunks("G", 0, 10)
    gw._require_consistent_chunks("G", 10, 10)
    with pytest.raises(RuntimeError, match="3 of 10"):
        gw._require_consistent_chunks("G", 3, 10)


def test_a_worker_that_cannot_load_its_gene_reports_it_and_does_not_return_an_empty_chunk(tmp_path):
    gene, result = gw._perms_worker_replay("GHOST", str(tmp_path), "no.nwk", None, "lg", str(tmp_path / "cache"), 0.1, ["b_1"], {"b_1": ([], [])}, "x")
    assert gene == "GHOST" and result is None


def test_a_worker_whose_asr_is_unavailable_reports_it(monkeypatch, tmp_path):
    monkeypatch.setattr(gw, "load_gene_context", lambda *a, **kw: None)     # the ASR could not be read or computed
    gene, result = gw._perms_worker_replay("PEPC", str(tmp_path), "t.nwk", None, "lg", str(tmp_path), 0.1, ["b_1"], {"b_1": ([], [])}, "x")
    assert gene == "PEPC" and result is None
