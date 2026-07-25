#!/usr/bin/env python3
# test_generation_golden.py — Golden-path generation tests vs. the real reference scripts.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Renders a ProjectConfig whose values mirror SBATCH_run_phenotypes_primates.sh /
run_phenotype_single_primates.sh and asserts semantic equivalence (SBATCH
directives, exported env vars, NF_FLAGS contents, CASE block arms) rather than
byte-identity — see implementation plan §9.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import subprocess

# ── Local ─────────────────────────────────────────────────────────────────────
from conftest import golden_project
from gui.generation.render import render_sbatch, render_single
from gui.generation.validate import validate


def _bash_syntax_ok(script: str) -> bool:
    result = subprocess.run(
        ["bash", "-n", "/dev/stdin"], input=script, text=True, capture_output=True
    )
    return result.returncode == 0


def test_golden_project_has_no_validation_errors():
    assert validate(golden_project()) == []


def test_sbatch_directives():
    out = render_sbatch(golden_project())
    assert "#SBATCH --job-name=phylophere" in out
    assert "#SBATCH --partition=haswell" in out
    assert "#SBATCH -t 144:00:00" in out
    assert "#SBATCH --mail-user=miguel.ramon@upf.edu" in out
    assert "#SBATCH --array=1-2%2" in out
    assert _bash_syntax_ok(out)


def test_sbatch_exports_module_toggles():
    out = render_sbatch(golden_project())
    assert "export RUN_REPORTING=true" in out
    assert "export RUN_CAAS=true" in out
    assert "export RUN_DISAMBIGUATION=true" in out
    assert "export RUN_ACCUMULATION=true" in out
    assert "export RUN_RER=false" in out
    assert "export RUN_FADE=false" in out
    assert "export RUN_VEP=true" in out
    assert "export RUN_SCORING=true" in out
    assert "export RUN_ENRICHMENT=true" in out
    assert "export RUN_POSENRICH=true" in out


def test_sbatch_phenotype_catalogue_case_block():
    out = render_sbatch(golden_project())
    assert 'CLASS=1; TRAIT="neoplasia_prevalence"' in out
    assert 'CLASS=1; TRAIT="malignant_prevalence"' in out
    assert 'export SCORING_RER_INPUT="' in out
    assert 'export SCORING_FADE_SUMMARY_TOP="' in out
    assert 'export SCORING_FADE_SUMMARY_BOTTOM="' in out
    assert _bash_syntax_ok(out)


def test_single_runner_class1_branch():
    out = render_single(golden_project())
    assert '--my_traits "$TRAIT_FILE"' in out
    assert '--branch_trait "$BRANCH_TRAIT"' in out
    assert '--prune_list "$PRUNE_LIST"' in out
    assert _bash_syntax_ok(out)


def test_single_runner_nf_flags_per_module():
    out = render_single(golden_project())
    assert "--contrast_selection" in out
    assert '--ct_tool "discovery,resample,bootstrap"' in out
    assert "--ct_disambiguation" in out
    assert "--ct_accumulation" in out
    assert "--vep" in out
    assert "--scoring" in out
    assert "--enrichment" in out
    assert "--posenrich" in out
    # RER/FADE disabled in the golden fixture, but their `if` blocks must still
    # exist (guarded by the RUN_RER/RUN_FADE env vars) so hand-edits stay possible.
    assert 'if [ "${RUN_RER:-false}" = true ]; then' in out
    assert 'if [ "${RUN_FADE:-false}" = true ]; then' in out
    assert _bash_syntax_ok(out)


def test_scoring_fade_site_fallbacks_are_exported():
    # Regression test: scoring_fade_site_top/bottom live on ScoringConfig (global,
    # non-per-row fallback) and must be exported by the SBATCH wrapper, unlike
    # scoring_fade_summary_top/bottom which are per-phenotype-row (see PhenotypeRow).
    project = golden_project()
    project.modules.scoring.scoring_fade_site_top = "/x/site_top.tsv"
    project.modules.scoring.scoring_fade_site_bottom = "/x/site_bottom.tsv"

    sbatch_out = render_sbatch(project)
    assert 'export SCORING_FADE_SITE_TOP="/x/site_top.tsv"' in sbatch_out
    assert 'export SCORING_FADE_SITE_BOTTOM="/x/site_bottom.tsv"' in sbatch_out

    single_out = render_single(project)
    assert '--scoring_fade_site_top "$SCORING_FADE_SITE_TOP"' in single_out
    assert '--scoring_fade_site_bottom "$SCORING_FADE_SITE_BOTTOM"' in single_out


def test_single_runner_nxf_home_collision_avoidance_preserved():
    out = render_single(golden_project())
    assert 'export NXF_HOME="$WORK_BASE/.nextflow/.nextflow_${TRAIT}"' in out
    assert 'ln -sfn "/homes/users/mramon/.nextflow/plugins" "$NXF_HOME/plugins"' in out
