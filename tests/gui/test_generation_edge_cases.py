#!/usr/bin/env python3
# test_generation_edge_cases.py — Disabled-module fallbacks, CLASS 2 rows, empty catalogue.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Standard library ──────────────────────────────────────────────────────────
import subprocess

# ── Local ─────────────────────────────────────────────────────────────────────
from conftest import golden_project
from gui.generation.render import render_batch, render_single
from gui.generation.validate import validate, validate_paths
from gui.models.runtime import PhenotypeRow


def _bash_syntax_ok(script: str) -> bool:
    result = subprocess.run(
        ["bash", "-n", "/dev/stdin"], input=script, text=True, capture_output=True
    )
    return result.returncode == 0


def test_caas_disabled_with_fallback_threads_through():
    project = golden_project()
    project.modules.caas.enabled = False
    project.precomputed.discovery_from = "/precomputed/discovery"
    project.precomputed.bootstrap_from = "/precomputed/bootstrap"

    assert validate(project) == []

    sbatch_out = render_batch(project)
    single_out = render_single(project)

    assert "export RUN_CAAS=false" in sbatch_out
    assert 'export DISCOVERY_FROM="/precomputed/discovery"' in sbatch_out
    assert 'export BOOTSTRAP_FROM="/precomputed/bootstrap"' in sbatch_out

    # The CAAS NF_FLAGS block stays behind a shell-runtime conditional (so a user
    # can still hand-flip RUN_CAAS after generation) — what must change is the
    # baked default the conditional falls back to, plus the fallback flags below.
    assert 'if [ "${RUN_CAAS:-false}" = true ]; then' in single_out
    assert '--discovery_from "$DISCOVERY_FROM"' in single_out
    assert '--bootstrap_from "$BOOTSTRAP_FROM"' in single_out
    assert _bash_syntax_ok(sbatch_out)
    assert _bash_syntax_ok(single_out)


def test_caas_disabled_without_fallback_fails_validation():
    project = golden_project()
    project.modules.caas.enabled = False
    # No discovery_from/resample_from/bootstrap_from supplied.

    errors = validate(project)
    assert any("CAAS is disabled" in e for e in errors)


def test_class2_row_renders_simple_trait_branch():
    project = golden_project()
    # RER/FADE are disabled in golden_project(); Scoring then requires a per-row
    # fallback for whichever row(s) exist, so supply one on this replacement row.
    project.runtime.phenotype_rows = [
        PhenotypeRow(
            trait_class=2,
            trait="diet_index",
            discrete_method="decile",
            scoring_rer_input="/data/rer_diet_index.tsv",
            scoring_fade_summary_top="/data/fade_top_diet_index.tsv",
            scoring_fade_summary_bottom="/data/fade_bottom_diet_index.tsv",
        ),
    ]
    project.runtime.simple_trait_file = "/data/diet_traitfile.csv"

    assert validate(project) == []

    sbatch_out = render_batch(project)
    single_out = render_single(project)

    assert 'CLASS=2; TRAIT="diet_index"' in sbatch_out
    assert 'DISCRETE="decile"' in sbatch_out
    assert '--my_traits "$SIMPLE_TRAIT_FILE"' in single_out
    assert '--discrete_method "$DISCRETE_METHOD"' in single_out
    # CLASS 2 branch must zero out the CLASS-1-only fields.
    assert '--n_trait ""' in single_out
    assert '--c_trait ""' in single_out
    assert '--secondary_trait ""' in single_out
    assert _bash_syntax_ok(sbatch_out)
    assert _bash_syntax_ok(single_out)


def test_empty_phenotype_catalogue_fails_validation():
    project = golden_project()
    project.runtime.phenotype_rows = []

    errors = validate(project)
    assert any("at least one row" in e for e in errors)


def test_module_disabled_without_fallback_flags_missing_downstream_input():
    project = golden_project()
    project.modules.rer.enabled = False
    project.modules.fade.enabled = False
    for row in project.runtime.phenotype_rows:
        row.scoring_rer_input = ""
        row.scoring_fade_summary_top = ""
        row.scoring_fade_summary_bottom = ""

    errors = validate(project)
    assert any("RERconverge is disabled" in e for e in errors)
    assert any("FADE is disabled" in e for e in errors)


def test_tower_disabled_bakes_false_default():
    project = golden_project()
    project.runtime.use_tower = False

    # -with-tower stays behind a shell-runtime conditional (hand-edit-friendly), so
    # what actually changes is the baked default the conditional falls back to.
    single_out = render_single(project)
    assert 'if [ "${RUN_TOWER:-false}" = true ]; then' in single_out
    assert "NF_FLAGS+=( -with-tower )" in single_out

    sbatch_out = render_batch(project)
    assert "export RUN_TOWER=false" in sbatch_out


def test_validate_paths_flags_nonexistent_paths():
    project = golden_project()  # fake /homes/users/mramon/... paths, don't exist here
    problems = validate_paths(project)
    assert problems  # every path is fake, so this must be non-empty
    assert any("repo directory" in p for p in problems)


def test_validate_paths_accepts_real_paths_and_ignores_empty():
    project = golden_project()
    project.general.repo_dir = "/tmp"
    project.general.nextflow_plugins_dir = "/tmp"
    project.runtime.ali_sp_names = ""  # empty -> must not be reported by validate_paths

    problems = validate_paths(project)
    assert not any("repo directory" in p for p in problems)
    assert not any("Nextflow plugins directory" in p for p in problems)
    assert not any("alignment species names" in p for p in problems)


def test_validate_paths_checks_phenotype_row_prune_files_relative_to_prune_dir():
    project = golden_project()
    project.runtime.prune_dir = "/tmp"
    project.runtime.phenotype_rows[0].prune = "definitely_missing_file.txt"

    problems = validate_paths(project)
    assert any("prune" in p and "definitely_missing_file.txt" in p for p in problems)
