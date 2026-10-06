#!/usr/bin/env python3
# project_io.py — Save/load a ProjectConfig to/from a JSON project file.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Project files: save and load a ProjectConfig as JSON (indent 2, trailing newline).

Thin wrapper around gui/models/serialization.py. It imports no PySide6: the file
dialogs belong to the widget layer, which calls these two functions with a path.

Imported by: gui/main.py, gui/widgets/main_window.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import json
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig
from gui.models.serialization import from_dict, to_dict


def save_project(path: Path, project: ProjectConfig) -> None:
    """Write `project` to `path` as JSON."""
    path.write_text(json.dumps(to_dict(project), indent=2) + "\n")


def load_project(path: Path) -> ProjectConfig:
    """Read the JSON project file at `path` (a schema mismatch raises ValueError)."""
    return from_dict(json.loads(path.read_text()))
