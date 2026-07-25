#!/usr/bin/env python3
# project_io.py — Save/load a ProjectConfig to/from a JSON project file.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Thin wrapper around gui/models/serialization.py. No PySide6 imports — QFileDialog
glue lives in the widget layer, which calls these two functions with a path.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import json
from pathlib import Path

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig
from gui.models.serialization import from_dict, to_dict


def save_project(path: Path, project: ProjectConfig) -> None:
    path.write_text(json.dumps(to_dict(project), indent=2) + "\n")


def load_project(path: Path) -> ProjectConfig:
    return from_dict(json.loads(path.read_text()))
