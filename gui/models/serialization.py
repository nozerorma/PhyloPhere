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
import logging
import typing
from typing import Any

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import SCHEMA_VERSION, ProjectConfig

logger = logging.getLogger(__name__)

# Batch-size parameters of the permulation null that ct_core_batch_size replaced; their values are not carried over.
_RETIRED_CAAS_FIELDS = ("ct_perm_replay_batch_size", "ct_disambig_perms_batch_size")

# Parameters of processes that no longer exist (the per-gene discovery task, the observed disambiguation chunking) or that
# nothing reads (the ASR mode: a gene is read from the ASR cache when it is there and computed when it is not).
_RETIRED_FIELDS = {"caas": ("ct_discovery_batch_size",), "disambiguation": ("ct_disambig_asr_mode", "ct_disambig_batch_size")}


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
    # The per-gene MAP directory moved from the VEP module (vep_map_dir) to post-processing (caas_map_dir).
    modules = data.get("modules")
    if isinstance(modules, dict):
        old = (modules.get("vep") or {}).pop("vep_map_dir", None)
        disambiguation = modules.get("disambiguation")
        if old and isinstance(disambiguation, dict) and not disambiguation.get("caas_map_dir"):
            disambiguation["caas_map_dir"] = old
        caas = modules.get("caas")
        if isinstance(caas, dict):
            if caas.pop("caas_b0_diagnostic", None) is not None:
                logger.warning("project parameter caas_b0_diagnostic was retired: the real labeling (b_0) is always replayed with the null")
            dropped = [k for k in _RETIRED_CAAS_FIELDS if caas.pop(k, None) is not None]
            if dropped:
                logger.warning("project parameters %s were replaced by ct_core_batch_size; the project uses its value (%s)",
                               ", ".join(dropped), caas.get("ct_core_batch_size", "default"))
        for module, names in _RETIRED_FIELDS.items():
            section = modules.get(module)
            gone = [k for k in names if isinstance(section, dict) and section.pop(k, None) is not None]
            if gone:
                logger.warning("project parameters %s were retired (module %s): no process reads them any more", ", ".join(gone), module)
    return data


def from_dict(data: dict[str, Any]) -> ProjectConfig:
    """Deserialize a project dict (already parsed from JSON) into a ProjectConfig."""
    data = migrate(data)
    return _dataclass_from_dict(ProjectConfig, data)
