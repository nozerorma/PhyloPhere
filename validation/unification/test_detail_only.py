"""`--detail-only` of the permulation null driver: pass A only, shards identical to the full run's.

The ASR replay workers are replaced by synthetic ones (module-level, so the pool pickles them by reference); the shard
writer, pass B0 and `_finalize_perm_scores` are the real code.
"""
import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows/CT_DISAMBIGUATION/local"))
import disambiguation_perms_main as dpm  # noqa: E402
from src.utils import gene_wrapper as gw  # noqa: E402

CYCLES = ["b_1", "b_2", "b_3"]
PASS_B_OUTPUTS = ["gene_cycle_scores.tsv", "perm_pos_cycle_caas.tsv.gz", "perm_pos_sample.tsv", "perm_pos_quantiles.tsv"]


def fake_replay(args):
    return args[0], [(c, []) for c in args[7]]


def fake_finalize(gene, pooled, n_cycles_total, postproc_filter, minlen, maxcaas, columns):
    rows = [{"Gene": gene, "cycle": c, "Position": 10 + i, "caap_group": "GS1",
             "asr_path_score": 0.25 * (i + 1) + 0.01 * len(gene), "n_detected": 1, "clust": 0, "side": "top"}
            for i, (c, _) in enumerate(pooled)]
    return gene, rows


@pytest.fixture
def lengths(tmp_path):
    f = tmp_path / "lengths.tsv"
    f.write_text("gene\tlength\nGENEA\t1000\nGENEB\t2000\n")
    return f


@pytest.fixture
def stubbed(monkeypatch):
    monkeypatch.setattr(gw, "build_cycle_inputs", lambda disc, resample, cycles: (list(CYCLES), {c: ([], []) for c in CYCLES}))
    monkeypatch.setattr(gw, "_perms_worker_replay_wrapper", fake_replay)
    monkeypatch.setattr(gw, "_perms_worker_finalize", fake_finalize)


def _run(out, lengths, **kw):
    # gene removal ON, so that pass B0 (removed_units.tsv) runs in the full run
    return gw.process_all_genes_perms(
        genes=["GENEA", "GENEB"], alignment_dir="x", tree_file="x", perm_discovery_file="x", resample_dir="x",
        taxid_mapping_path=None, asr_model="lg", asr_cache_dir="x", posterior_threshold=0.0, workers=1,
        output_dir=out, postproc_filter=True, gene_lengths_file=str(lengths), gene_filter_mode="extreme", **kw)


def _shards(out):
    return {p.name: gzip.open(p, "rt").read() for p in sorted((out / "perm_pos_detail").glob("*.tsv.gz"))}


def test_detail_only_stops_after_pass_a(stubbed, lengths, tmp_path):
    out = _run(tmp_path / "o", lengths, detail_only=True)
    assert set(_shards(out)) == {"GENEA.tsv.gz", "GENEB.tsv.gz"}
    assert (out / "perm_pos_detail.manifest.tsv").exists()
    for name in PASS_B_OUTPUTS + ["removed_units.tsv"]:
        assert not (out / name).exists(), name


def test_detail_only_shards_equal_those_of_the_full_run(stubbed, lengths, tmp_path):
    full = _run(tmp_path / "full", lengths)
    only = _run(tmp_path / "only", lengths, detail_only=True)
    assert all((full / n).exists() for n in PASS_B_OUTPUTS + ["removed_units.tsv"])  # the default still runs passes B0 and B
    assert _shards(only) == _shards(full) and _shards(full)


def test_the_cli_flag_reaches_the_null_and_the_b0_runs(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(dpm, "process_all_genes_perms", lambda **kw: calls.append(kw) or kw["output_dir"])
    monkeypatch.setattr(dpm, "_genes_from_export", lambda p: (["GENEA"], {"GENEA": 1}))
    monkeypatch.setattr(dpm, "_read_resample_labelings", lambda d: {"b_0": 0, "b_1": 0})
    base = ["x", "--alignment-dir", "a", "--tree", "t", "--perm-discovery", "d", "--resample-dir", "r",
            "--output-dir", str(tmp_path / "o"), "--asr-cache-dir", "c"]
    for extra, expected in (([], False), (["--detail-only"], True)):
        calls.clear()
        monkeypatch.setattr(sys, "argv", base + extra)
        dpm.main()
        assert len(calls) == 2  # null + b_0
        assert [bool(c.get("detail_only")) for c in calls] == [expected, expected]
