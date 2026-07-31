#!/usr/bin/env python3
# render.py — Jinja2 rendering: ProjectConfig -> the two generated shell scripts.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure functions, no PySide6 import — reusable headless (tests, a future CLI) without
pulling in Qt. See gui/generation/context.py for the shared env-var naming contract
between the two templates.
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

    single_runner_filename overrides the auto-computed name of the companion
    per-phenotype script this one dispatches to (SINGLE_RUNNER=...) — must match
    whatever filename that script is actually saved as (see
    MainWindow.generate_scripts, which overrides RuntimeConfig.script_base_name),
    or the generated batch script looks for a file that doesn't exist.
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
