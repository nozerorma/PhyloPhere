#!/usr/bin/env python3
# render.py — Jinja2 rendering: ProjectConfig → the batch and single-phenotype shell scripts.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure functions, no PySide6 import, so scripts can be rendered without Qt. The batch
script (sbatch_array.sh.j2) exports the run configuration and dispatches one
phenotype row per call to the single-phenotype script (run_single.sh.j2). The
environment uses StrictUndefined: a template variable missing from the context
raises instead of rendering empty. See gui/generation/context.py for the context
and the variable-name contract between the two templates.

Imported by: gui/widgets/main_window.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from functools import lru_cache

# ── Third-party ───────────────────────────────────────────────────────────────
import jinja2

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.generation.context import build_context
from gui.models.project import ProjectConfig


@lru_cache(maxsize=1)
def _environment() -> jinja2.Environment:
    # Built once; trim_blocks and lstrip_blocks keep the {% %} control lines out of
    # the generated shell, and keep_trailing_newline preserves the final newline.
    return jinja2.Environment(
        loader=jinja2.PackageLoader("gui.generation", "templates"),
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
        undefined=jinja2.StrictUndefined,
    )


def render_batch(
    project: ProjectConfig,
    postproc_mode: str | None = None,
    reuse_exploratory: bool = False,
    single_runner_filename: str | None = None,
) -> str:
    """Render the batch-runner script for the given project.

    postproc_mode ("exploratory" or "filter") selects the post-processing pass the
    pair of scripts runs; reuse_exploratory makes the "filter" pass read the output
    of the exploratory pass instead of recomputing it.

    single_runner_filename overrides the derived name of the per-phenotype script
    this one dispatches to (SINGLE_RUNNER). It must match the filename under which
    that script is saved (MainWindow.generate_scripts applies
    RuntimeConfig.script_base_name to both), or the batch script points at a file
    that does not exist.
    """
    ctx = build_context(project)
    ctx["reuse_exploratory"] = reuse_exploratory
    ctx["postproc_mode"] = postproc_mode or ""
    ctx["single_runner_filename"] = ""
    if postproc_mode:
        ctx["postproc_mode"] = postproc_mode
        ctx["single_runner_filename"] = f"run_phenotype_single_{postproc_mode if postproc_mode != 'filter' else 'complete'}.sh"
        if "modules" in ctx and "disambiguation" in ctx["modules"]:
            ctx["modules"]["disambiguation"]["caas_postproc_mode"] = postproc_mode
    if single_runner_filename is not None:
        ctx["single_runner_filename"] = single_runner_filename
    return _environment().get_template("sbatch_array.sh.j2").render(**ctx)


def render_single(
    project: ProjectConfig, postproc_mode: str | None = None, reuse_exploratory: bool = False
) -> str:
    """Render the single-phenotype runner script for the given project."""
    ctx = build_context(project)
    ctx["reuse_exploratory"] = reuse_exploratory
    if postproc_mode:
        ctx["postproc_mode"] = postproc_mode
        if "modules" in ctx and "disambiguation" in ctx["modules"]:
            ctx["modules"]["disambiguation"]["caas_postproc_mode"] = postproc_mode
    return _environment().get_template("run_single.sh.j2").render(**ctx)
