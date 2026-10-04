"""Stage 2.0 wiring: the null's gene universe, parameter defaults, the MAP directory and the cluster-file selection.

The tests run Nextflow itself (`-preview` for the DAG, small scripts run with `-main-script` for defaults, helper
functions and the CT_FILTER process) and are skipped when `nextflow` is not on PATH. PHYLOPHERE_ROOT points them at
another checkout (used to show that a test fails on the tree before a change).
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
MINI = HERE  # the mini scripts live next to this file, in the tree under test
needs_nextflow = pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")


# ── helpers ──────────────────────────────────────────────────────────────────

def _dag_edges(mmd_text):
    """Labelled-node edges of a Nextflow mermaid DAG; anonymous operator nodes are collapsed."""
    label = {}
    for m in re.finditer(r'^\s*(v\d+)\s*[\[\(\{]+(?:"([^"]*)")?', mmd_text, re.M):
        label[m.group(1)] = (m.group(2) or "").strip()
    succ = {}
    for m in re.finditer(r'^\s*(v\d+)\s*-->(?:\|[^|]*\|)?\s*(v\d+)', mmd_text, re.M):
        succ.setdefault(m.group(1), set()).add(m.group(2))
    named = {n for n, l in label.items() if l}
    edges = set()
    for a in named:
        stack, seen = list(succ.get(a, ())), set()
        while stack:
            b = stack.pop()
            if b in seen:
                continue
            seen.add(b)
            if b in named:
                edges.add((label[a], label[b]))
            else:
                stack.extend(succ.get(b, ()))
    return edges


def _preview(tmp_path, **overrides):
    """DAG edges of main.nf with a minimal CT + null + post-processing configuration."""
    r = _run_preview(tmp_path, **overrides)
    assert (tmp_path / "dag.mmd").exists(), r.stdout[-600:] + r.stderr[-600:]
    return _dag_edges((tmp_path / "dag.mmd").read_text())


def _run_preview(tmp_path, **overrides):
    """The finished `nextflow -preview` process of main.nf under the same configuration (the DAG is at tmp_path/dag.mmd)."""
    (tmp_path / "ali").mkdir()
    shutil.copy(ROOT / "validation/unification/golden/pepc_c4_complete/PEPC.fasta", tmp_path / "ali/PEPC.Homo_sapiens.fa")
    (tmp_path / "traits").mkdir()
    (tmp_path / "traits/traitfile_H1.tab").write_text("a\t1\nb\t0\n")
    (tmp_path / "tree.nwk").write_text("(a:1,b:1);\n")
    (tmp_path / "values.tab").write_text("a\t1\n")
    (tmp_path / "genes.tsv").write_text("gene\tlength\nPEPC\t1000\n")
    params = {"outdir": str(tmp_path / "out"), "reporting": False, "alignment": str(tmp_path / "ali"),
              "tree": str(tmp_path / "tree.nwk"), "my_traits": str(tmp_path / "values.tab"),
              "caas_config": str(tmp_path / "traits"), "gene_ensembl_file": str(tmp_path / "genes.tsv"),
              "ct_tool": "discovery,resample", "ct_disambiguation": True, "ct_postproc": True,
              "caas_permulation_enrichment": True, "ct_core_batch_size": "20", "caas_full_perms": "10",
              "seed": "1998"}
    params.update(overrides)
    (tmp_path / "params.json").write_text(json.dumps(params))
    r = subprocess.run(["nextflow", "run", str(ROOT / "main.nf"), "-preview", "-profile", "local", "-params-file",
                        str(tmp_path / "params.json"), "-with-dag", str(tmp_path / "dag.mmd")],
                       cwd=tmp_path, capture_output=True, text=True, timeout=240)
    return r


def _mini(tmp_path, script, *params):
    """Run a mini script as the main.nf of a project whose folders link the tree under test (baseDir resolves there)."""
    proj = tmp_path / "proj"
    proj.mkdir(exist_ok=True)
    for name in ("subworkflows", "workflows", "conf", "lib", "bin", "nextflow.config"):
        if (ROOT / name).exists() and not (proj / name).exists():
            (proj / name).symlink_to(ROOT / name)
    shutil.copy(MINI / script, proj / "main.nf")
    # the pipeline sets PYTHONNOUSERSITE=1 for container isolation; a local task without a container needs the user site
    (tmp_path / "test.config").write_text("env { PYTHONNOUSERSITE = '' }\n")
    return subprocess.run(["nextflow", "run", str(proj / "main.nf"), "-profile", "local", "-c", str(tmp_path / "test.config"), *params],
                          cwd=tmp_path, capture_output=True, text=True, timeout=240)


# ── the null's gene universe ─────────────────────────────────────────────────

@needs_nextflow
def test_the_null_universe_is_the_cleaned_background_when_the_null_is_batched(tmp_path):
    edges = _preview(tmp_path)
    assert ("CAAS_BACKGROUND_CLEANUP", "CAAS_CORE_MERGE") in edges


@needs_nextflow
def test_the_null_universe_is_the_cleaned_background_when_the_null_is_not_batched(tmp_path):
    edges = _preview(tmp_path, ct_core_batch_size="1")
    assert ("CAAS_BACKGROUND_CLEANUP", "CAAS_CORE_MERGE") in edges


_NULL_ONLY = dict(ct_tool="", enrichment=True, caas_permulation_enrichment=True, ct_disambiguation=False, ct_postproc=False)
_OLD_NULL_PROCESSES = ("PERM_REPLAY", "PERM_REPLAY_BATCHED", "CAAS_PERMS_DISAMBIGUATE", "CAAS_PERMS_DISAMBIGUATE_BATCHED",
                       "CAAS_PERMS_AGGREGATE", "CAAS_PERMS_MERGE_DETAIL", "CAAS_PERMS_REBUILD")


@needs_nextflow
def test_the_live_null_replays_and_disambiguates_in_one_process_family(tmp_path):
    edges = _preview(tmp_path)
    nodes = {n for e in edges for n in e}
    assert not nodes & set(_OLD_NULL_PROCESSES)
    assert {("SUBSET_RESAMPLE_PERMS", "CAAS_CORE_BATCHED"), ("Channel.fromList", "CAAS_CORE_BATCHED"),
            ("CAAS_CORE_BATCHED", "CAAS_CORE_MERGE"), ("CAAS_CORE_MERGE", "CT_ACCUMULATION_RANDOMIZE")} <= edges


_RETIRED_OBSERVED = ("DISCOVERY", "DISCOVERY_BATCHED", "CONCAT_DISCOVERY", "CONCAT_BACKGROUND", "CT_DISAMBIGUATION_RUN",
                     "CT_DISAMBIGUATION_SPLIT_GENES", "CT_DISAMBIGUATION_RUN_BATCHED", "CT_DISAMBIGUATION_MERGE",
                     "CT_DISAMBIGUATION_PLOTS", "CAAS_PERMULATION")


@needs_nextflow
def test_the_live_observed_labeling_is_the_b0_slice_of_the_core_and_feeds_everything_downstream(tmp_path):
    edges = _preview(tmp_path)
    assert not {n for e in edges for n in e} & set(_RETIRED_OBSERVED)
    # the core's b_0 slice -> the contract files -> the report, the post-processing, the universe, ASR robustness
    assert {("CAAS_CORE_BATCHED", "CAAS_CORE_OBSERVED"), ("CAAS_CORE_OBSERVED", "CAAS_META_CAAS_REPORT"),
            ("CAAS_CORE_OBSERVED", "CAAS_PREPARE_POSTPROC_INPUT"), ("CAAS_CORE_OBSERVED", "CAAS_BACKGROUND_CLEANUP"),
            ("CAAS_CORE_OBSERVED", "ASR_ROBUSTNESS_REPORT")} <= edges
    # the null merge waits for the cleaned universe, which needs the observed master: the observed step cannot wait for the merge
    assert ("CAAS_BACKGROUND_CLEANUP", "CAAS_CORE_MERGE") in edges and ("CAAS_CORE_MERGE", "CAAS_CORE_OBSERVED") not in edges
    assert ("CAAS_CORE_MERGE", "CAAS_PREPARE_POSTPROC_INPUT") not in edges


@needs_nextflow
def test_a_discovery_tab_that_already_exists_is_scored_by_caas_observed_and_the_core_writes_no_observed_files(tmp_path):
    (tmp_path / "disc").mkdir()
    (tmp_path / "disc/discovery.tab").write_text("gene\tmode\tcaap_group\ttrait\tposition\tcaas\tamino_encoded\tpattern\nPEPC\tCAAP\tUS\tH1\t3\tA/B\tA/B\t1\n")
    (tmp_path / "disc/background_genes.output").write_text("PEPC\n")
    edges = _preview(tmp_path, ct_tool="resample", discovery_from=str(tmp_path / "disc/discovery.tab"),
                     background_input=str(tmp_path / "disc/background_genes.output"))
    nodes = {n for e in edges for n in e}
    assert "CAAS_OBSERVED" in nodes and "CAAS_CORE_OBSERVED" not in nodes and not nodes & set(_RETIRED_OBSERVED)
    assert ("CAAS_OBSERVED", "CAAS_PREPARE_POSTPROC_INPUT") in edges and ("CAAS_OBSERVED", "ASR_ROBUSTNESS_REPORT") in edges


@needs_nextflow
def test_discovery_without_the_permuted_labelings_says_what_is_missing(tmp_path):
    r = _run_preview(tmp_path, ct_tool="discovery", caas_permulation_enrichment=False)
    assert "needs the permuted labelings" in r.stdout + r.stderr


@needs_nextflow
def test_perm_discovery_exports_that_already_exist_go_to_the_same_process(tmp_path):
    pub = tmp_path / "out/caas_permulation"
    (pub / "perm_disc").mkdir(parents=True)
    (pub / "resample_perms.tab").write_text("b_0~H1\ta\tb\nb_1~H1\ta\tb\n")
    (pub / "perm_disc/G.perm_replay.discovery.output").write_text("x\n")
    edges = _preview(tmp_path, **_NULL_ONLY)
    assert ("Channel.fromPath", "CAAS_CORE_BATCHED") in edges and ("Channel.fromList", "CAAS_CORE_BATCHED") not in edges
    assert not {n for e in edges for n in e} & ({"SUBSET_RESAMPLE_PERMS"} | set(_OLD_NULL_PROCESSES))


@needs_nextflow
@pytest.mark.parametrize("first_row", ["b_1~H1\ta\tb\n", "b_1\ta\tb\n"])
def test_perm_discovery_exports_made_without_the_real_labeling_are_rejected(tmp_path, first_row):
    pub = tmp_path / "out/caas_permulation"
    (pub / "perm_disc").mkdir(parents=True)
    (pub / "resample_perms.tab").write_text(first_row)
    (pub / "perm_disc/G.perm_replay.discovery.output").write_text("x\n")
    r = _run_preview(tmp_path, **_NULL_ONLY)
    assert r.returncode != 0 and "without the real labeling (b_0)" in r.stdout + r.stderr


@needs_nextflow
def test_the_standalone_null_subsets_the_resample_and_replays_the_alignments(tmp_path):
    (tmp_path / "out/caastools").mkdir(parents=True)
    (tmp_path / "out/caastools/resample.tab").write_text("b_1\ta\tb\n")
    edges = _preview(tmp_path, **_NULL_ONLY)
    assert {("caas_config", "SUBSET_RESAMPLE_PERMS"), ("SUBSET_RESAMPLE_PERMS", "CAAS_CORE_BATCHED"),
            ("Channel.fromList", "CAAS_CORE_BATCHED"), ("caas_config", "CAAS_CORE_BATCHED")} <= edges


def test_main_sends_each_source_of_the_null_to_its_own_input_of_caas_core():
    """Source-level check: the DAG collapses which input of CAAS_CORE_BATCHED a channel feeds."""
    text = (ROOT / "main.nf").read_text()
    for line in ("perm_align_ch = ct_results.caas_align_tuple", "perm_reuse_ch = Channel.fromPath(precomp_disc_files).collect()",
                 "perm_align_ch = align_tuple_standalone"):
        assert line in text, line
    call = text[text.index("core = CAAS_CORE("):].split(")")[0]
    assert [a.strip() for a in call.split("(")[1].split(",")][:3] == ["perm_align_ch", "perm_reuse_ch", "perm_cfg_ch"]


def test_main_feeds_the_observed_step_with_the_b0_slices_and_the_merge_with_the_shards():
    """Source-level check: the DAG cannot tell which output of CAAS_CORE reaches which consumer."""
    text = (ROOT / "main.nf").read_text()
    assert "CAAS_CORE_OBSERVED(core.b0_observed.collect()," in text
    assert "CAAS_CORE_MERGE(core.pos_detail.collect()," in text
    assert text.index("CAAS_CORE_OBSERVED(") < text.index("CT_POSTPROC(") < text.index("CAAS_CORE_MERGE(")  # the universe needs the observed master


def _publish_dirs(path, process):
    """[(directory, pattern)] of the publishDir lines of one process definition."""
    text = (ROOT / path).read_text()
    body = text[text.index(f"process {process} {{"):]
    body = body[:body.index("    input:")]
    return [(m.group(1), m.group(2)) for m in re.finditer(r'publishDir path: "\$\{params\.outdir\}([^"]*)"[^\n]*?pattern: \'([^\']*)\'', body)]


def test_the_observed_processes_publish_their_output_directories_by_name():
    """A process that declares a directory as output publishes it by its name: `dir/**` matches no declared output and
    publishes nothing while the run succeeds (the files in outdir are checked by test_core_batched.py)."""
    expected = {("/caastools", "{discovery.tab,background.output,background_genes.output}"), ("/meta_caas", "meta_caas"),
                ("", "ct_disambiguation")}
    assert set(_publish_dirs("subworkflows/CT/caas_permulation.nf", "CAAS_CORE_OBSERVED")) == expected
    assert set(_publish_dirs("subworkflows/CT_DISAMBIGUATION/ct_observed.nf", "CAAS_OBSERVED")) == {("/meta_caas", "meta_caas"), ("", "ct_disambiguation")}


@needs_nextflow
def test_accumulation_with_the_permulation_null_stops_early_when_the_null_has_no_cycle(tmp_path):
    r = _run_preview(tmp_path, caas_full_perms="0", ct_accumulation=True, accumulation_randomization_type="permulation")
    text = r.stdout + r.stderr
    assert r.returncode != 0 and "caas_full_perms 0" in text and "permulation" in text, text[-800:]


@needs_nextflow
def test_accumulation_with_another_randomization_runs_without_permuted_labelings(tmp_path):
    edges = _preview(tmp_path, caas_full_perms="0", ct_accumulation=True, accumulation_randomization_type="cons_decile")
    assert ("CAAS_CORE_OBSERVED", "CAAS_CORE_MERGE") not in edges and any(dst == "CAAS_CORE_MERGE" for _, dst in edges)


@needs_nextflow
def test_accumulation_with_the_permulation_null_runs_when_there_are_permuted_labelings(tmp_path):
    edges = _preview(tmp_path, caas_full_perms="10", ct_accumulation=True, accumulation_randomization_type="permulation")
    assert any(dst == "CAAS_CORE_MERGE" for _, dst in edges)


def test_the_standalone_scoring_route_rebuilds_the_null_with_the_merge_process():
    """Source-level check (the scoring route needs inputs `-preview` does not accept); the process itself is run by
    test_core_batched.py on one shard directory, the shape this route passes."""
    text = (ROOT / "workflows/scoring.nf").read_text()
    assert "include { CAAS_CORE_MERGE }" in text
    assert "_rebuild = CAAS_CORE_MERGE(Channel.value(_detail), resolved_background," in text
    assert "CAAS_PERMS_REBUILD" not in text


def test_the_null_merge_treats_any_no_prefixed_universe_as_absent():
    text = (ROOT / "subworkflows/CT/caas_permulation.nf").read_text()
    assert "universe.name != 'NO_FILE'" not in text
    assert len(re.findall(r"universe\.name\.startsWith\('NO_'\)", text)) == 1


# ── defaults and the MAP parameter ───────────────────────────────────────────

def _params_line(r):
    m = re.search(r"PARAMS filter_maxcaas=(\S+) caas_map_dir=\[(.*)\]", r.stdout)
    assert m, r.stdout[-400:] + r.stderr[-400:]
    return m.group(1), m.group(2)


@needs_nextflow
def test_filter_maxcaas_has_a_default_and_caas_map_dir_is_empty_by_default(tmp_path):
    maxcaas, map_dir = _params_line(_mini(tmp_path, "mini_params.nf"))
    assert float(maxcaas) == 0.7 and map_dir == ""


@needs_nextflow
def test_caas_map_dir_reads_an_older_vep_map_dir_and_the_new_name_wins(tmp_path):
    assert _params_line(_mini(tmp_path, "mini_params.nf", "--vep_map_dir", "/old"))[1] == "/old"
    assert _params_line(_mini(tmp_path, "mini_params.nf", "--vep_map_dir", "/old", "--caas_map_dir", "/new"))[1] == "/new"


# ── the null's post-processing arguments ─────────────────────────────────────

def _functions(tmp_path, *params):
    r = _mini(tmp_path, "mini_functions.nf", *params)
    return {k: v for k, v in re.findall(r"^([A-Z_0-9]+)=(.*)$", r.stdout, re.M)}, r


@needs_nextflow
def test_the_null_receives_the_map_directory_only_when_it_is_set(tmp_path):
    off, r = _functions(tmp_path)
    assert "--postproc-filter" in off["ARGS"] and "--clust-maxcaas 0.7" in off["ARGS"], r.stdout[-500:] + r.stderr[-300:]
    assert "--train-map-dir" not in off["ARGS"]
    on, _ = _functions(tmp_path, "--caas_map_dir", "/maps")
    assert "--train-map-dir /maps" in on["ARGS"]


@needs_nextflow
def test_the_null_postproc_arguments_are_empty_when_the_filters_are_off_or_the_annotation_is_a_sentinel(tmp_path):
    got, _ = _functions(tmp_path)
    assert got["ARGS_SENTINEL"] == "[]" and got["ON"] == "true"
    got, _ = _functions(tmp_path, "--caas_perms_postproc", "false")
    assert got["ARGS"] == "[]" and got["ON"] == "false"


# ── exploratory grid and the selected cluster file ───────────────────────────

@needs_nextflow
def test_the_exploratory_grid_holds_the_selected_pair(tmp_path):
    got, r = _functions(tmp_path)
    inside = got["GRID_IN"]                      # the selected pair (3, 0.7) is in the 4 x 3 grid
    assert inside.count("[") == 1 + 12 and "[3, 0.7]" in inside, inside + r.stderr[-300:]
    outside = got["GRID_OUT"]                    # the selected pair (5, 0.65) is not in the 3 x 2 grid and is added
    assert outside.count("], [") == 7 - 1 and outside.endswith("[5, 0.65]]"), outside


@needs_nextflow
def test_the_cluster_file_suffix_matches_the_name_the_filter_script_writes(tmp_path):
    got, _ = _functions(tmp_path)
    (tmp_path / "in.tsv").write_text("Gene\tPosition\tcaap_group\n" + "".join(f"G\t{p}\tUS\n" for p in (1, 2, 3, 40)))
    script = ROOT / "subworkflows/CT_POSTPROC/local/src/filter_caas_clusters-param.py"
    for minlen, maxcaas, key in ((3, 0.7, "SUFFIX_70"), (2, 0.29, "SUFFIX_29")):
        subprocess.run([sys.executable, str(script), "-i", str(tmp_path / "in.tsv"), "-l", str(minlen), "-c", str(maxcaas)],
                       check=True, capture_output=True)
        written = [p.name for p in tmp_path.glob("in.filtered.*.tsv") if p.name.endswith(got[key])]
        assert written, (got[key], [p.name for p in tmp_path.glob("in.*")])


# ── the real CT_FILTER process with the MAP directory ────────────────────────

def _run_ct_filter(tmp_path, with_map):
    """Positions that CT_FILTER discards for gene A, whose trimmed positions 9, 10, 11 sit at untrimmed columns 10, 14, 20."""
    tmp_path.mkdir()
    (tmp_path / "in.tsv").write_text("Gene\tPosition\tcaap_group\n" + "".join(f"A\t{p}\tUS\n" for p in (9, 10, 11, 40)))
    maps = tmp_path / "maps"
    maps.mkdir()
    selected = list(range(1, 10)) + [10, 14, 20] + list(range(21, 60))
    sel = {c: i + 1 for i, c in enumerate(selected)}
    with open(maps / "A.Lemur_catta.map.tsv", "w") as fh:
        fh.write("ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\n")
        for c in range(1, 60):
            fh.write(f"{c}\t{'selected' if c in sel else 'removed'}\t{sel.get(c, 'NA')}\t{sel.get(c, 'NA')}\n")
    extra = ["--caas_map_dir", str(maps)] if with_map else []
    r = _mini(tmp_path, "mini_ct_filter.nf", "--mini_input", str(tmp_path / "in.tsv"), "--outdir", str(tmp_path / "out"), *extra)
    out = tmp_path / "out/postproc/filter_selected/in.filtered.minlen3.maxcaas70.tsv"
    assert out.exists(), r.stdout[-600:] + r.stderr[-600:]
    rows = [l.split("\t") for l in out.read_text().splitlines()[1:]]
    return {int(pos) for gene, pos, group, flag in rows if flag == "Discarded"}


@needs_nextflow
def test_ct_filter_without_a_map_directory_measures_trains_in_trimmed_positions(tmp_path):
    assert _run_ct_filter(tmp_path / "off", with_map=False) == {9, 10, 11}


@needs_nextflow
def test_ct_filter_measures_trains_in_untrimmed_columns_when_the_map_directory_is_set(tmp_path):
    assert _run_ct_filter(tmp_path / "on", with_map=True) == set()
