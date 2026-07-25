#!/usr/bin/env python3
# serialization.py — Generic dataclass <-> dict (<-> JSON) round-trip for ProjectConfig.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Field-declaration order gives stable, human-diffable JSON key ordering "for free"
via dataclasses.asdict(). Reconstruction (`from_dict`) is type-hint-driven so nested
dataclasses (GeneralConfig, RuntimeConfig, ModulesConfig, ResourcesConfig,
list[PhenotypeRow], and each of the 8 module configs) round-trip without per-class
boilerplate.

`migrate()` is a no-op dispatch stub today (schema_version is always 1) but exists
so a future schema change doesn't require a breaking rewrite of project_io.py.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import dataclasses
import typing
from typing import Any

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import SCHEMA_VERSION, ProjectConfig


def to_dict(project: ProjectConfig) -> dict[str, Any]:
    """Serialize a ProjectConfig to a plain, JSON-ready, key-order-stable dict."""
    return dataclasses.asdict(project)


def _reconstruct(field_type: Any, value: Any) -> Any:
    origin = typing.get_origin(field_type)
    if origin is list:
        (item_type,) = typing.get_args(field_type)
        return [_reconstruct(item_type, item) for item in value]
    if dataclasses.is_dataclass(field_type):
        return _dataclass_from_dict(field_type, value)
    return value


def _dataclass_from_dict(cls: type, data: dict[str, Any]) -> Any:
    hints = typing.get_type_hints(cls)
    kwargs = {}
    for f in dataclasses.fields(cls):
        if f.name not in data:
            continue  # let the dataclass's own default/default_factory apply
        kwargs[f.name] = _reconstruct(hints[f.name], data[f.name])
    return cls(**kwargs)


def migrate(data: dict[str, Any]) -> dict[str, Any]:
    """Upgrade an older on-disk project dict to the current schema, if needed.

    No-op today (only schema_version 1 exists). Future migrations should branch on
    data.get("schema_version") and mutate/return a new dict at the current version.
    """
    version = data.get("schema_version", SCHEMA_VERSION)
    if version != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported project schema_version={version!r}; "
            f"this build only supports version {SCHEMA_VERSION}."
        )
    return data


def from_dict(data: dict[str, Any]) -> ProjectConfig:
    """Deserialize a project dict (already parsed from JSON) into a ProjectConfig."""
    data = migrate(data)
    return _dataclass_from_dict(ProjectConfig, data)
