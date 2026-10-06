#!/usr/bin/env python3
# resource_defaults.py — Reads conf/resources.config into ProcessResourceOverride rows
# for the "load conf defaults" button of the Resources tab.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Resource defaults: parse the per-process cpus and memory of conf/resources.config into
override rows.

The module imports no PySide6, so it can be tested without a display like
gui/generation/.

The parser is regex-based, not a Groovy parser. It relies on the regular shape of
conf/resources.config: one `withName:` or `withLabel:` selector per block, `cpus = { N }`
and `memory = { ... }` as single-line assignments, and no nested braces inside a selector
block. Brace depth is tracked, so a block closes when the depth returns to where it opened
and not at the first bare "}" line. A block whose selector is an alternation
(`withName: 'A|B'`) is not a row and is skipped, as is a block with neither cpus nor memory
(errorStrategy-only labels).

conf/resources.config is the only source of per-process defaults. The override table starts
empty, and a row there is a deliberate deviation from it; this parser only fills the table
with the current defaults as a starting point to edit.

`memory` keeps the whole closure body, not only the leading `N.GB`: several entries are
`memory = { N.GB * task.attempt }`, which gives a retry more memory on each attempt. The rows
are rendered verbatim into the `-c` override config (run_single.sh.j2), so cutting the
`* task.attempt` would pin the process at its first-attempt memory and defeat the retry.

Imported by: gui/widgets/resource_table/widget.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import re
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.resources import ProcessResourceOverride

CONF_DIR = Path(__file__).resolve().parent.parent / "conf"

DEFAULTS_FILE = CONF_DIR / "resources.config"

_SELECTOR_RE = re.compile(r"with(Name|Label):?\s*'?([A-Za-z0-9_]+)'?\s*\{")
_CPUS_RE = re.compile(r"cpus\s*=\s*\{\s*(\d+)\s*\}")
# Whole closure body (e.g. "16.GB * task.attempt"), not only the leading number and unit.
_MEM_RE = re.compile(r"memory\s*=\s*\{\s*([^}]+?)\s*\}")


def load_defaults(path: Path = DEFAULTS_FILE) -> list[ProcessResourceOverride]:
    """Parse `path` (default conf/resources.config) into a new list of override rows."""
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
