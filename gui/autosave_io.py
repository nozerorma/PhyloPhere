#!/usr/bin/env python3
# autosave_io.py — Persists the live (possibly-unsaved) project state across restarts.
# PhyloPhere | gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Distinct from project_io.py's explicit Save/Open: this tracks whatever is
currently in memory — dirty or not — so quitting the GUI (or it crashing)
never silently loses in-progress edits the user never got around to saving.
Written on every project mutation and read back once at startup
(see gui/widgets/main_window.py's _write_autosave/__init__).

Best-effort by design: a write failure (e.g. no disk space) or a corrupt/
unreadable file on read must never block the UI or crash startup — worst case
is losing the autosave, not the app.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import json
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtCore import QStandardPaths

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig
from gui.models.serialization import from_dict, to_dict

AUTOSAVE_FILENAME = "autosave.json"


def _autosave_path() -> Path:
    app_data_dir = Path(
        QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    )
    app_data_dir.mkdir(parents=True, exist_ok=True)
    return app_data_dir / AUTOSAVE_FILENAME


def write_autosave(project: ProjectConfig, project_path: Path | None, dirty: bool) -> None:
    payload = {
        "project_path": str(project_path) if project_path else None,
        "dirty": dirty,
        "project": to_dict(project),
    }
    try:
        _autosave_path().write_text(json.dumps(payload))
    except OSError:
        pass


def read_autosave() -> tuple[ProjectConfig, Path | None, bool] | None:
    """Returns (project, project_path, dirty) restored from the last autosave, or
    None if there isn't one / it can't be parsed (treated as "start fresh")."""
    path = _autosave_path()
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
        project = from_dict(payload["project"])
        project_path = Path(payload["project_path"]) if payload.get("project_path") else None
        dirty = bool(payload.get("dirty", False))
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return project, project_path, dirty
