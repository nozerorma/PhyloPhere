"""CAAS_EVIDENCE on PEPC with precomputed ASR: Nextflow runs the real process, the reference is explain_positions.py run by hand.

The evidence of the N best positions is published under <outdir>/scoring/evidence/. The process is off unless
`caas_evidence_top_n` is above 0, and a run that cannot give it the observed discovery.tab says so before it starts.
"""
import gzip
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402

GOLD = tw.ROOT / "validation/unification/golden/pepc_c4_complete"
LOCAL = tw.ROOT / "subworkflows/CT_DISAMBIGUATION/local"


@pytest.fixture(scope="module")
def inp(tmp_path_factory):
    d = tmp_path_factory.mktemp("evidence")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    (d / "ali").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", d / "ali/PEPC.fasta")
    with gzip.open(GOLD / "discovery.tab.gz", "rt") as src, open(d / "discovery.tab", "w") as dst:
        shutil.copyfileobj(src, dst)
    return d


def _params(inp, out, top, cache="asr_cache"):
    """Command line of the mini script; an empty cache dir goes in a params file (`--x ""` reaches Nextflow as `true`)."""
    i = inp / "observed_inputs"
    p = ["--mini_discovery", str(inp / "discovery.tab"), "--mini_scores", str(GOLD / "position_scores.tsv"),
         "--mini_design", str(i / "traitfiles"), "--mini_tree", str(i / "pruned_tree_file.nwk"), "--outdir", str(out),
         "--alignment", str(inp / "ali"), "--tax_id", str(i / "taxid.tsv"), "--gene_ensembl_file", str(i / "gene_ensembl.tsv"),
         "--ct_disambig_posterior_threshold", "0.1", "--ct_disambig_asr_model", "lg", "--caas_evidence_top_n", str(top)]
    if cache:
        return p + ["--ct_disambig_asr_cache_dir", str(i / cache)]
    pf = out.parent / "empty_cache.json"
    pf.parent.mkdir(parents=True, exist_ok=True)
    pf.write_text('{"ct_disambig_asr_cache_dir": ""}')
    return p + ["-params-file", str(pf)]


@tw.needs_nextflow
def test_the_published_evidence_is_what_the_script_gives_by_hand(inp, tmp_path):
    r = tw._mini(tmp_path, "mini_caas_evidence.nf", *_params(inp, tmp_path / "out", 3))
    listing = tmp_path / "out/evidence_paths.txt"
    assert listing.exists(), r.stdout[-1500:] + r.stderr[-1500:]
    published = tmp_path / "out/scoring/evidence"
    assert sorted(f.name for f in published.iterdir()) == ["evidence_top3.tsv", "top_positions.tsv"]
    i = inp / "observed_inputs"
    subprocess.run([sys.executable, str(LOCAL / "explain_positions.py"), "--alignment-dir", str(inp / "ali"), "--tree", str(i / "pruned_tree_file.nwk"),
                    "--discovery", str(inp / "discovery.tab"), "--position-scores", str(GOLD / "position_scores.tsv"), "--top", "3",
                    "--design", str(i / "traitfiles"), "--output-dir", str(tmp_path / "hand"), "--asr-model", "lg", "--posterior-threshold", "0.1",
                    "--workers", "2", "--asr-cache-dir", str(i / "asr_cache"), "--taxid-mapping", str(i / "taxid.tsv"),
                    "--ensembl-genes-file", str(i / "gene_ensembl.tsv")], check=True, capture_output=True)
    for name in ("evidence_top3.tsv", "top_positions.tsv"):
        assert (published / name).read_bytes() == (tmp_path / "hand" / name).read_bytes()
    assert (published / "evidence_top3.tsv").read_text().count("\n") > 12


@tw.needs_nextflow
def test_an_empty_asr_cache_dir_stops_the_process(inp, tmp_path):
    r = tw._mini(tmp_path, "mini_caas_evidence.nf", *_params(inp, tmp_path / "out", 3, cache=""))
    # Nextflow echoes the whole script in its error report, so the message is read from the task's own stderr
    errs = [f.read_text() for f in (tmp_path / "work").glob("*/*/.command.err")]
    assert any("ct_disambig_asr_cache_dir must be set" in e for e in errs), r.stdout[-800:] + r.stderr[-800:]
    assert not (tmp_path / "out/scoring/evidence").exists()


@tw.needs_nextflow
def test_ct_observed_emits_the_design_and_tree_it_scored_with(inp, tmp_path):
    """CAAS_EVIDENCE explains with the design and tree of the scoring it follows: CT_OBSERVED must hand them on."""
    p = _params(inp, tmp_path / "out", 0)
    i = inp / "observed_inputs"
    r = tw._mini(tmp_path, "mini_ct_observed_emits.nf", *p[:p.index("--outdir")], "--outdir", str(tmp_path / "out"), *p[p.index("--alignment"):])
    emits = tmp_path / "out/emits.txt"
    assert emits.exists(), r.stdout[-1200:] + r.stderr[-1200:]
    assert emits.read_text().splitlines() == [f"design {i / 'traitfiles'}", f"tree {i / 'pruned_tree_file.nwk'}"]


# ── the wiring in main.nf ────────────────────────────────────────────────────

_SCORED = dict(traitname="t", scoring=True, vep=False, enrichment=False)


def _scoring_inputs(tmp_path):
    f = tmp_path / "filtered_discovery.tsv"
    f.write_text("Gene\tPosition\n")
    return f


@tw.needs_nextflow
def test_the_process_is_off_by_default(tmp_path):
    edges = tw._preview(tmp_path, **_SCORED)
    assert ("CAAS_CORE_MERGE", "SCORING_COMPUTE") in edges
    assert not any("CAAS_EVIDENCE" in e for e in edges)


@tw.needs_nextflow
def test_a_live_run_feeds_it_the_scores_and_the_observed_discovery_of_the_core(tmp_path):
    edges = tw._preview(tmp_path, caas_evidence_top_n=2, **_SCORED)
    assert ("SCORING_COMPUTE", "CAAS_EVIDENCE") in edges
    assert ("CAAS_CORE_OBSERVED", "CAAS_EVIDENCE") in edges


@tw.needs_nextflow
def test_a_run_that_scores_a_given_discovery_feeds_it_that_discovery(tmp_path):
    disc = tmp_path / "given" / "caastools" / "discovery.tab"
    disc.parent.mkdir(parents=True)
    disc.write_text("gene\tposition\n")
    (disc.parent / "background_genes.output").write_text("PEPC\n")
    edges = tw._preview(tmp_path, ct_tool="", discovery_from=str(disc), background_input=str(disc.parent / "background_genes.output"),
                        caas_evidence_top_n=2, **_SCORED)
    assert ("SCORING_COMPUTE", "CAAS_EVIDENCE") in edges
    assert ("CAAS_CORE_OBSERVED", "CAAS_EVIDENCE") not in edges   # the discovery is the given file, not the core's


@tw.needs_nextflow
def test_a_run_without_the_observed_discovery_is_refused_before_it_starts(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    only_scoring = dict(ct_tool="", ct_disambiguation=False, ct_postproc=False, scoring_postproc_input=str(_scoring_inputs(tmp_path)), **_SCORED)
    ok = tw._run_preview(tmp_path / "a", caas_evidence_top_n=0, **only_scoring)
    assert ok.returncode == 0, ok.stdout[-800:] + ok.stderr[-800:]
    r = tw._run_preview(tmp_path / "b", caas_evidence_top_n=2, **only_scoring)
    assert r.returncode != 0 and "caas_evidence_top_n" in r.stdout + r.stderr


@tw.needs_nextflow
def test_the_evidence_needs_scoring(tmp_path):
    r = tw._run_preview(tmp_path, caas_evidence_top_n=2, traitname="t", vep=False, enrichment=False)
    assert r.returncode != 0 and "caas_evidence_top_n" in r.stdout + r.stderr
