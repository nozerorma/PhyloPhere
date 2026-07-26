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


def render_batch(project: ProjectConfig) -> str:
    """Render the batch-runner script for the given project: a SLURM array-job
    wrapper when runtime.runtime_type == "slurm", or a plain sequential bash loop
    over every phenotype row when "local" (see sbatch_array.sh.j2's own internal
    {% if profile == 'slurm' %} branching — one template, two execution modes,
    sharing the same per-phenotype dispatch logic either way)."""
    ctx = build_context(project)
    return _environment().get_template("sbatch_array.sh.j2").render(**ctx)


def render_single(project: ProjectConfig) -> str:
    """Render the single-phenotype runner script for the given project."""
    ctx = build_context(project)
    return _environment().get_template("run_single.sh.j2").render(**ctx)
