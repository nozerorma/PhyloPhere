#!/usr/bin/env python3
# project.py — ProjectConfig: the single source-of-truth container for the whole GUI.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Every tab reads/writes a slice of this object; gui/generation/ renders it to the two
shell scripts; gui/project_io.py (de)serializes it to a JSON project file.

`schema_version` and `flavor` exist from day one even though there is currently only
one schema version and one flavor ("primates") — this keeps a future TOGA-style
env-var-only pathway (explicitly deferred, see implementation plan) an additive
change rather than a breaking migration.
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
