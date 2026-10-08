#!/usr/bin/env python3
# context.py — ProjectConfig → Jinja2 render context for the two shell templates.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure function, no PySide6 import. Flattens a ProjectConfig into the dict both
templates render from: the config sections as plain dicts, the module toggles
resolved to the strings "true"/"false" (run_defaults), the CT and RER tool strings,
the runtime profile and the SLURM array size and spec.

The shell variable names themselves are written literally in sbatch_array.sh.j2
(which exports them) and run_single.sh.j2 (which reads them, with run_defaults as
fallback), so a rename must be made in both templates.

Imported by: gui/generation/render.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import dataclasses
from typing import Any

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig


def _bool_str(value: bool) -> str:
    # Shell and JSON boolean literal, as the templates write it.
    return "true" if value else "false"


def _ct_tool_string(caas) -> str:
    # Comma-separated --ct_tool value ("discovery,resample") from the CAAS tool toggles.
    parts = []
    if caas.ct_tool_discovery:
        parts.append("discovery")
    if caas.ct_tool_resample:
        parts.append("resample")
    return ",".join(parts)


def _rer_tool_string(rer) -> str:
    # Comma-separated --rer_tool value from the RERconverge tool toggles.
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
    """Build the full Jinja2 render context for both templates.

    ctx keeps the config sections (general, runtime, modules, resources, precomputed)
    as plain dicts and adds run_defaults, ct_tool_string, rer_tool_string, profile,
    array_size and array_spec.
    """
    ctx = dataclasses.asdict(project)

    caas = project.modules.caas
    disambig = project.modules.disambiguation
    scoring = project.modules.scoring
    enrichment = project.modules.enrichment

    ctx["ct_tool_string"] = _ct_tool_string(caas)
    ctx["rer_tool_string"] = _rer_tool_string(project.modules.rer)

    ctx["profile"] = project.runtime.runtime_type  # "local" or "slurm"

    pc = project.precomputed
    # A "use precomputed X" box on the Precomputed Run tab wins over that module's own
    # enabled toggle. The tab's cascade (PrecomputedTab._toggle_module) keeps the two
    # mutually exclusive, but it only fires on an interactive checkbox click: a project
    # loaded from a hand-edited JSON or a template, or one whose tabs were out of sync
    # at generation time, can carry enabled=true and use_x=true together. Rendering
    # that combination would run the stage live and also feed it a precomputed result,
    # so the precedence is enforced here, once, instead of being assumed upstream.
    # Post-processing runs with Disambiguation, so it follows disambig.enabled.
    caas_enabled = caas.enabled and not (pc.use_discovery or pc.use_resample or pc.use_ct)
    disambiguation_enabled = disambig.enabled and not pc.use_disambiguation
    ct_postproc_enabled = disambig.enabled and not pc.use_postproc
    accumulation_enabled = project.modules.accumulation.enabled and not pc.use_accumulation
    rer_enabled = project.modules.rer.enabled and not pc.use_rer
    fade_enabled = project.modules.fade.enabled and not pc.use_fade
    vep_enabled = project.modules.vep.enabled and not pc.use_vep

    ctx["run_defaults"] = {
        "reporting": _bool_str(project.general.reporting),
        "tower": _bool_str(project.runtime.use_tower),
        "toy_mode": _bool_str(project.runtime.toy_mode),
        "caas": _bool_str(caas_enabled),
        "caas_permulation_enrichment": _bool_str(caas.caas_permulation_enrichment),
        "caas_perms_postproc": _bool_str(caas.caas_perms_postproc),
        "disambiguation": _bool_str(disambiguation_enabled),
        "ct_postproc": _bool_str(ct_postproc_enabled),
        "run_postproc_exploratory": _bool_str(getattr(disambig, 'run_postproc_exploratory', True)),
        "run_postproc_filter": _bool_str(getattr(disambig, 'run_postproc_filter', True)),
        "accumulation": _bool_str(accumulation_enabled),
        "rer": _bool_str(rer_enabled),
        "fade": _bool_str(fade_enabled),
        "vep": _bool_str(vep_enabled),
        "scoring": _bool_str(scoring.enabled),
        "enrichment": _bool_str(enrichment.enabled),
        "posenrich": _bool_str(enrichment.posenrich_enabled),
        "scoring_ami": _bool_str(getattr(enrichment, 'scoring_ami', enrichment.scoring_string)),
        "scoring_string": _bool_str(enrichment.scoring_string),
        "publish_domino_intermediates": _bool_str(enrichment.publish_domino_intermediates),
        "auto_generate_ensembl": _bool_str(scoring.auto_generate_ensembl),
        "fcs_enabled": _bool_str(enrichment.fcs_enabled),
        "auto_fetch_gmt": _bool_str(enrichment.auto_fetch_gmt),
        "auto_fetch_eggnog": _bool_str(enrichment.auto_fetch_eggnog),
        "posenrich_domains": _bool_str(enrichment.posenrich_domains),
        "posenrich_ucr": _bool_str(enrichment.posenrich_ucr),
        "posenrich_eggnog": _bool_str(enrichment.posenrich_eggnog),
        "miss_pair": _bool_str(caas.miss_pair),
        "caap_mode": _bool_str(caas.caap_mode),
        "resample_use_n": _bool_str(getattr(caas, 'resample_use_n', True)),
        "multi_hypothesis": _bool_str(getattr(caas, 'multi_hypothesis', True)),
        "perm_match_pss": _bool_str(getattr(caas, 'perm_match_pss', True)),
        "publish_intermediates": _bool_str(caas.publish_intermediates),
        "asr_diagnostics": _bool_str(disambig.asr_diagnostics),
    }

    # One array task per phenotype row; the optional concurrency cap becomes "%N".
    n_rows = len(project.runtime.phenotype_rows)
    cap = project.runtime.sbatch_array_concurrency.strip()
    ctx["array_size"] = n_rows
    ctx["array_spec"] = f"1-{n_rows}" + (f"%{cap}" if cap else "")

    return ctx
