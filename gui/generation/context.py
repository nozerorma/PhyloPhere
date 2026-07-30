#!/usr/bin/env python3
# context.py — ProjectConfig -> Jinja2 render-context dict.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure function, no PySide6 import. Single source of truth for the shell env-var
names shared between sbatch_array.sh.j2 (which `export`s them) and run_single.sh.j2
(which reads them) — keep the two templates' variable names in sync with this file.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import dataclasses
from typing import Any

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig


def _bool_str(value: bool) -> str:
    return "true" if value else "false"


def _ct_tool_string(caas) -> str:
    parts = []
    if caas.ct_tool_discovery:
        parts.append("discovery")
    if caas.ct_tool_resample:
        parts.append("resample")
    if caas.ct_tool_bootstrap:
        parts.append("bootstrap")
    return ",".join(parts)


def _rer_tool_string(rer) -> str:
    parts = []
    if rer.rer_tool_build_trait:
        parts.append("build_trait")
    if rer.rer_tool_build_tree:
        parts.append("build_tree")
    if rer.rer_tool_build_matrix:
        parts.append("build_matrix")
    if rer.rer_tool_continuous:
        parts.append("continuous")
    return ",".join(parts)


def build_context(project: ProjectConfig) -> dict[str, Any]:
    """Build the full Jinja2 render context for both templates."""
    ctx = dataclasses.asdict(project)  # general/runtime/modules/resources as plain dicts

    caas = project.modules.caas
    disambig = project.modules.disambiguation
    scoring = project.modules.scoring
    enrichment = project.modules.enrichment

    ctx["ct_tool_string"] = _ct_tool_string(caas)
    ctx["rer_tool_string"] = _rer_tool_string(project.modules.rer)

    ctx["profile"] = project.runtime.runtime_type  # "local" | "slurm"

    ctx["run_defaults"] = {
        "reporting": _bool_str(project.general.reporting),
        "tower": _bool_str(project.runtime.use_tower),
        "toy_mode": _bool_str(project.runtime.toy_mode),
        "caas": _bool_str(caas.enabled),
        "caas_permulation_enrichment": _bool_str(caas.caas_permulation_enrichment),
        "disambiguation": _bool_str(disambig.enabled),
        "ct_postproc": _bool_str(disambig.ct_postproc),
        "run_postproc_exploratory": _bool_str(getattr(disambig, 'run_postproc_exploratory', True)),
        "run_postproc_filter": _bool_str(getattr(disambig, 'run_postproc_filter', True)),
        "accumulation": _bool_str(project.modules.accumulation.enabled),
        "rer": _bool_str(project.modules.rer.enabled),
        "fade": _bool_str(project.modules.fade.enabled),
        "vep": _bool_str(project.modules.vep.enabled),
        "scoring": _bool_str(scoring.enabled),
        "scoring_stress": _bool_str(scoring.scoring_stress),
        "enrichment": _bool_str(enrichment.enabled),
        "posenrich": _bool_str(enrichment.posenrich_enabled),
        "scoring_ami": _bool_str(getattr(enrichment, 'scoring_ami', enrichment.scoring_string)),
        "scoring_string": _bool_str(enrichment.scoring_string),
        "rer_permulation_enrichment": _bool_str(enrichment.rer_permulation_enrichment),
        "miss_pair": _bool_str(caas.miss_pair),
        "caap_mode": _bool_str(caas.caap_mode),
        "include_b0": _bool_str(caas.include_b0),
        "publish_intermediates": _bool_str(caas.publish_intermediates),
        "export_groups": _bool_str(caas.export_groups),
        "export_perm_discovery": _bool_str(caas.export_perm_discovery),
        "asr_robustness": _bool_str(disambig.asr_robustness),
    }

    n_rows = len(project.runtime.phenotype_rows)
    cap = project.runtime.sbatch_array_concurrency.strip()
    ctx["array_size"] = n_rows
    ctx["array_spec"] = f"1-{n_rows}" + (f"%{cap}" if cap else "")

    return ctx
