#!/usr/bin/env python3
# serialization.py — Dataclass to dict conversion of ProjectConfig for its JSON file.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Serialization: ProjectConfig to a JSON-ready dict (`to_dict`) and back (`from_dict`).

dataclasses.asdict() keeps the field-declaration order, so the JSON keys come out in a
stable order that diffs cleanly. Reconstruction follows the type hints, which rebuilds
the nested dataclasses (GeneralConfig, RuntimeConfig, ModulesConfig, ResourcesConfig,
list[PhenotypeRow] and each module config) without per-class code. Keys absent from the
dict take the dataclass default.

`migrate()` runs first in `from_dict`: it rejects an unsupported schema_version, renames
vep.vep_map_dir to disambiguation.caas_map_dir, and drops the parameters that no
process reads, with a logged warning.

Imported by: gui/project_io.py, gui/autosave_io.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import dataclasses
import logging
import typing
from typing import Any

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import SCHEMA_VERSION, ProjectConfig

logger = logging.getLogger(__name__)

# Batch-size parameters of the permulation null that ct_core_batch_size covers; the value of ct_core_batch_size is used instead.
_RETIRED_CAAS_FIELDS = ("ct_perm_replay_batch_size", "ct_disambig_perms_batch_size")

# cleanup_policy of older project files -> (retired intermediate flag, remove_work)
_CLEANUP_POLICY_OPTIONS = {"none": (False, False), "prune_intermediates": (True, False), "remove_work": (True, True)}

_RETIRED_PROCESS_SELECTORS = {
    "process_perm_replay", "process_perm_replay_batched", "CAAS_PERMS_DISAMBIGUATE",
    "CAAS_PERMS_DISAMBIGUATE_BATCHED", "CAAS_PERMS_AGGREGATE", "CAAS_PERMS_MERGE_DETAIL", "CAAS_PERMS_REBUILD",
}

# Parameters that no process reads: those of tasks the pipeline does not have (the per-gene discovery task, the chunking of
# the observed disambiguation) and the ASR mode (a gene is read from the ASR cache when it is there and computed when not).
_RETIRED_FIELDS = {
    "caas": ("ct_discovery_batch_size", "export_groups", "export_perm_discovery"),
    "disambiguation": ("ct_disambig_asr_mode", "ct_disambig_batch_size"),
    "scoring": ("scoring_weight_caas", "scoring_weight_rer", "scoring_weight_fade", "scoring_rer_direction", "scoring_stress", "scoring_stress_top_n",
                "scoring_stress_rank_metric", "scoring_ami", "scoring_string", "scoring_compare_fdr", "scoring_compare_top_n"),
    "enrichment": ("posenrich_eggnog", "cosmic_db", "fcs_fdr_wilcoxon", "fcs_fdr_lachenbruch"),
    "accumulation": ("accumulation_n_randomizations", "accumulation_randomization_type", "accumulation_entropy_dir"),
}


def to_dict(project: ProjectConfig) -> dict[str, Any]:
    """Convert a ProjectConfig to a plain dict of JSON types, with stable key order."""
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
            continue  # the dataclass default or default_factory applies
        kwargs[f.name] = _reconstruct(hints[f.name], data[f.name])
    return cls(**kwargs)


def migrate(data: dict[str, Any]) -> dict[str, Any]:
    """Bring a project dict read from disk to the current layout.

    Raises ValueError when schema_version is not SCHEMA_VERSION. A missing
    schema_version is taken as current. The dict is modified in place and returned.
    """
    version = data.get("schema_version", SCHEMA_VERSION)
    if version != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported project schema_version={version!r}; "
            f"this build only supports version {SCHEMA_VERSION}."
        )
    scoring = data.setdefault("modules", {}).setdefault("scoring", {})
    rer = data.get("modules", {}).get("rer", {})
    for name in ("scoring_gene_top_pct", "scoring_position_top_pct"):
        if name in scoring:
            scoring.pop(name)
            logger.warning("Retired %s: percentiles are replaced by empirical p-values; its value is not a p-value threshold", name)
    if "rer_pval_threshold" in rer:
        scoring.setdefault("evidence_p_threshold", rer.pop("rer_pval_threshold"))
    if "rer_pval_column" in rer:
        rer.pop("rer_pval_column")
        logger.warning("Retired rer_pval_column: reports show nominal and BH empirical evidence together")
    # Migrate cleanup policies while preserving the explicit work-removal setting.
    runtime = data.get("runtime")
    if isinstance(runtime, dict) and "cleanup_policy" in runtime:
        policy = runtime.pop("cleanup_policy")
        if policy not in _CLEANUP_POLICY_OPTIONS:
            logger.warning("project parameter cleanup_policy=%r is not a known value; remove_work is left off", policy)
        _retired_intermediates, remove_work = _CLEANUP_POLICY_OPTIONS.get(policy, (False, False))
        runtime.setdefault("remove_work", remove_work)
    if isinstance(runtime, dict) and "remove_intermediates" in runtime:
        runtime.pop("remove_intermediates")
        logger.warning("project parameter remove_intermediates was retired: task hygiene is unconditional; published null blocks are preserved")
    resources = data.get("resources")
    if isinstance(resources, dict) and isinstance(resources.get("process_overrides"), list):
        rows = resources["process_overrides"]
        resources["process_overrides"] = [row for row in rows if row.get("selector") not in _RETIRED_PROCESS_SELECTORS]
        if len(rows) != len(resources["process_overrides"]):
            logger.warning("resource overrides for retired permulation processes were dropped")
    # Older project files hold the per-gene MAP directory as vep.vep_map_dir; it is
    # copied to disambiguation.caas_map_dir unless that is already set.
    modules = data.get("modules")
    if isinstance(modules, dict):
        old = (modules.get("vep") or {}).pop("vep_map_dir", None)
        disambiguation = modules.get("disambiguation")
        if old and isinstance(disambiguation, dict) and not disambiguation.get("caas_map_dir"):
            disambiguation["caas_map_dir"] = old
        caas = modules.get("caas")
        if isinstance(caas, dict):
            old_stat = caas.get("gene_stat", "max")
            if old_stat != "max":
                logger.warning("Retired CAAS gene_stat=%s: max is now the only statistic; legacy gene scores/null matrices must be rebuilt", old_stat)
                caas["gene_stat"] = "max"
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
    """Build a ProjectConfig from a project dict already parsed from JSON."""
    data = migrate(data)
    return _dataclass_from_dict(ProjectConfig, data)
