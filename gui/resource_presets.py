#!/usr/bin/env python3
# resource_presets.py — Parses conf/resources.config.{slurm,local,local_lowspec}
# into ProcessResourceOverride rows for the Resources tab's preset buttons.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
No PySide6 import — hermetic and unit-testable like gui/generation/, even though
it lives one level up (it's shared by both the Resources tab widget and, in
principle, script generation, not generation-specific).

Regex-based, not a real Groovy parser: relies on conf/resources.config's existing
regularity (one `withName:`/`withLabel:` selector per block, `cpus = { N }` and
`memory = { N.GB }` as single-line assignments, no nested braces inside a
selector block) — see gui/generation/templates' own comment about keeping these
files' shape in sync. Brace-depth tracked explicitly rather than assumed, so a
selector block closes when depth returns to where it opened, not on the first
bare "}" line.

`memory` captures the FULL closure body, not just the leading `N.GB` — several
presets use `memory = { N.GB * task.attempt }` so a retry (see errorStrategy in
conf/resources.config) gets more memory on each attempt. An earlier version of
this regex captured only the numeric+unit prefix, silently truncating
`* task.attempt` off every scaled entry; since these rows get rendered verbatim
into a `-c`-loaded override config that is generated once per run and never
revisited, that truncation permanently pinned every overridden process at its
attempt-1 memory for the life of the project, defeating retry-with-more-memory
without any error or warning (a real production incident, not a hypothetical:
`POSENRICH_BUILD_GMT` OOM'd on attempt 1 and kept OOMing identically on retry).
"""

# ── Standard library ──────────────────────────────────────────────────────────
import re
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.resources import ProcessResourceOverride

CONF_DIR = Path(__file__).resolve().parent.parent / "conf"

PRESET_FILES = {
    "slurm": CONF_DIR / "resources.config.slurm",
    "local": CONF_DIR / "resources.config.local",
    "local_lowspec": CONF_DIR / "resources.config.local_lowspec",
}

_SELECTOR_RE = re.compile(r"with(Name|Label):?\s*'?([A-Za-z0-9_]+)'?\s*\{")
_CPUS_RE = re.compile(r"cpus\s*=\s*\{\s*(\d+)\s*\}")
# Full closure body (e.g. "16.GB * task.attempt"), not just the leading
# numeric+unit — see module docstring.
_MEM_RE = re.compile(r"memory\s*=\s*\{\s*([^}]+?)\s*\}")


def load_preset(name: str) -> list[ProcessResourceOverride]:
    """Parses one of the three bundled presets into a fresh list of override rows."""
    path = PRESET_FILES[name]
    rows: list[ProcessResourceOverride] = []
    current: ProcessResourceOverride | None = None
    depth = 0
    open_depth = 0
    for line in path.read_text().splitlines():
        sel = _SELECTOR_RE.search(line)
        if sel and current is None:
            current = ProcessResourceOverride(
                selector_type="withName" if sel.group(1) == "Name" else "withLabel",
                selector=sel.group(2),
            )
            open_depth = depth
        if current is not None:
            m = _CPUS_RE.search(line)
            if m:
                current.cpus = m.group(1)
            m = _MEM_RE.search(line)
            if m:
                current.memory = m.group(1)
        depth += line.count("{") - line.count("}")
        if current is not None and depth == open_depth:
            if current.cpus or current.memory:  # skip errorStrategy-only labels
                rows.append(current)
            current = None
    return rows
