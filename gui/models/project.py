#!/usr/bin/env python3
# project.py — ProjectConfig: the container that holds the whole state of the GUI.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
ProjectConfig: every tab reads and writes one slice of this object, gui/generation/
renders it into the two shell scripts, and gui/project_io.py (de)serializes it to a
JSON project file.

`schema_version` identifies the layout of the JSON file (checked by
serialization.migrate). `flavor` names the dataset family of the project; "primates"
is the only value.

Imported by: gui/project_io.py, gui/autosave_io.py, gui/generation/, gui/widgets/main_window.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field
from typing import Literal

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.general import GeneralConfig
from gui.models.modules import ModulesConfig
from gui.models.precomputed import PrecomputedConfig
from gui.models.resources import ResourcesConfig
from gui.models.runtime import RuntimeConfig

SCHEMA_VERSION = 1


@dataclass(kw_only=True)
class ProjectConfig:
    schema_version: int = SCHEMA_VERSION
    flavor: Literal["primates"] = "primates"
    general: GeneralConfig = field(default_factory=GeneralConfig)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    modules: ModulesConfig = field(default_factory=ModulesConfig)
    precomputed: PrecomputedConfig = field(default_factory=PrecomputedConfig)
    resources: ResourcesConfig = field(default_factory=ResourcesConfig)
