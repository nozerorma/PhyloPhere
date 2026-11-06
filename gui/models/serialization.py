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
vep.vep_map_dir to disambiguation.caas_map_dir, replaces disambiguation.gene_filter_mode by the two gene-removal switches, and drops the parameters that no
process reads, with a logged warning.

Imported by: gui/project_io.py, gui/autosave_io.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
import dataclasses
import logging
import re
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
    "disambiguation": ("ct_disambig_asr_mode", "ct_disambig_batch_size", "exploratory_without_null"),
    "scoring": ("scoring_weight_caas", "scoring_weight_rer", "scoring_weight_fade", "scoring_rer_direction", "scoring_stress", "scoring_stress_top_n",
                "scoring_stress_rank_metric"),
    "enrichment": ("posenrich_eggnog", "cosmic_db", "fcs_fdr_wilcoxon", "fcs_fdr_lachenbruch"),
    "posenrich": ("posenrich_ranking", "eggnog_taxid", "auto_fetch_eggnog", "egg_members_file", "egg_annotations_file"),
    "fade": ("fade_min_genes_for_heatmap",),
}

# Fields of the POSENRICH module that older projects hold in the enrichment section.
_POSENRICH_FIELDS = (
    "accumulation_caas_input", "accumulation_background_input",
    "posenrich_domains", "domain_variability_file", "domain_ref_species", "pfam_cache_dir",
    "posenrich_ucr", "ucr_positions_file", "fubar_sites_file",
    "posenrich_radios", "posenrich_custom_file",
)
# Memory parameters of POSENRICH that became per-process resource overrides (Resources tab): parameter -> process selector.
_POSENRICH_MEMORY_OVERRIDES = {
    "posenrich_within_gene_memory": "POSENRICH_WITHIN_GENE",
    "position_registry_memory": "POSITION_REGISTRY",
    "position_module_memory": "POSENRICH_CHARACTERIZE|POSENRICH_PROXIMITY",
}
# Parameters of the position-set enrichment test, which no process reads.
_RETIRED_POSENRICH_FIELDS = ("accumulation_fdr", "posenrich_min_size", "posenrich_max_size", "posenrich_batch_size", "posenrich_enabled",
                             "posenrich_n_permutations", "posenrich_seed", "posenrich_aggregation", "position_annotation_memory",
                             "posenrich_ranking", "eggnog_taxid", "auto_fetch_eggnog", "egg_members_file", "egg_annotations_file")


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


def _memory_params_to_overrides(data: dict[str, Any], *sections: Any) -> None:
    """POSENRICH memory parameters of older projects (kept in the posenrich section, or in the enrichment one of the oldest) become
    rows of the per-process resource overrides (a row already there wins)."""
    resources = None
    for name, selector in _POSENRICH_MEMORY_OVERRIDES.items():
        values = [str(section.pop(name, "") or "").strip() for section in sections if isinstance(section, dict)]
        value = next((v for v in values if v), "")
        if not value:
            continue
        match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([KMGT]B)", value, re.I)
        expression = f"{match.group(1)}.{match.group(2).upper()}" if match else value
        resources = data.setdefault("resources", {}) if resources is None else resources
        rows = resources.setdefault("process_overrides", [])
        if any(row.get("selector") == selector for row in rows):
            logger.warning("project parameter %s was replaced by a process override; the existing override of %s is kept", name, selector)
            continue
        rows.append({"selector_type": "withName", "selector": selector, "cpus": "", "memory": expression})
        logger.warning("project parameter %s became a process override of %s (memory %s) in the Resources tab", name, selector, expression)


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
    if "evidence_fdr_threshold" in scoring:
        scoring.pop("evidence_fdr_threshold")
        logger.warning("Retired evidence_fdr_threshold: BH evidence is classified at evidence_fdr_levels (default 0.05,0.1,0.15)")
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
    if isinstance(runtime, dict) and runtime.pop("clade_taxid", None) is not None:
        logger.warning("project parameter clade_taxid was retired: no process reads it")
    if isinstance(runtime, dict) and "remove_intermediates" in runtime:
        runtime.pop("remove_intermediates")
        logger.warning("project parameter remove_intermediates was retired: task hygiene is unconditional; published null blocks are preserved")
    posenrich = data["modules"].get("posenrich")
    if isinstance(posenrich, dict) and "posenrich_holdout" in posenrich:
        posenrich.pop("posenrich_holdout")
        logger.warning("project parameter posenrich_holdout was retired: no process reads it")
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
        # POSENRICH has its own section. Older projects hold its fields in the enrichment section (and the Accumulation
        # FDR in a module of its own); a value already in the posenrich section wins.
        posenrich = modules.setdefault("posenrich", {})
        old_acc = modules.pop("accumulation", None)
        retired_fdr = posenrich.pop("accumulation_fdr", None) is not None
        if retired_fdr or (isinstance(old_acc, dict) and "accumulation_fdr" in old_acc):
            logger.warning("Retired accumulation_fdr: the POSENRICH report classifies q at evidence_fdr_levels (default 0.05,0.1,0.15)")
        enrichment = modules.get("enrichment")
        if isinstance(enrichment, dict):
            if enrichment.pop("posenrich_modular", False):
                posenrich.setdefault("enabled", True)
            for name in _POSENRICH_FIELDS:
                if name in enrichment:
                    posenrich.setdefault(name, enrichment.pop(name))
            gone = [k for k in _RETIRED_POSENRICH_FIELDS if enrichment.pop(k, None) is not None]
            if gone:
                logger.warning("project parameters %s were retired (module posenrich): no process reads them any more", ", ".join(gone))
        _memory_params_to_overrides(data, posenrich, modules.get("enrichment"))
        scoring_section = modules.get("scoring")
        if isinstance(scoring_section, dict):
            evidence_file = scoring_section.pop("caas_evidence_discovery", "")
            if evidence_file:
                flag = f"--caas_evidence_discovery {evidence_file}"
                scoring_section["extra_flags"] = (str(scoring_section.get("extra_flags", "")).strip() + " " + flag).strip()
                logger.warning("project parameter caas_evidence_discovery left the GUI (the generated scripts set it when they reuse tables); "
                               "its value was added to the extra flags of the Scoring tab")
        precomputed = data.get("precomputed")
        if isinstance(precomputed, dict):
            precomputed.pop("use_accumulation", None)
        old = (modules.get("vep") or {}).pop("vep_map_dir", None)
        disambiguation = modules.get("disambiguation")
        if old and isinstance(disambiguation, dict) and not disambiguation.get("caas_map_dir"):
            disambiguation["caas_map_dir"] = old
        # gene_filter_mode (none|stratified) became two switches, one per rule.
        if isinstance(disambiguation, dict) and "gene_filter_mode" in disambiguation:
            mode = disambiguation.pop("gene_filter_mode")
            if mode == "none":
                disambiguation["remove_extreme_genes"] = disambiguation["remove_dubious_genes"] = False
                logger.warning("Retired gene_filter_mode=none: remove_extreme_genes and remove_dubious_genes are now both off")
            else:
                if mode != "stratified":
                    logger.warning("Retired gene_filter_mode=%r is not a known value: both gene-removal switches are left on", mode)
                disambiguation["remove_extreme_genes"] = disambiguation["remove_dubious_genes"] = True
        caas = modules.get("caas")
        if isinstance(caas, dict):
            old_stat = caas.get("gene_stat", "max")
            if old_stat != "max":
                logger.warning("Retired CAAS gene_stat=%s: max is now the only statistic; legacy gene scores/null matrices must be rebuilt", old_stat)
                caas["gene_stat"] = "max"
            if caas.pop("strict_conservation_output", None) is not None:
                logger.warning("project parameter strict_conservation_output was retired: the strict tables are written exactly when 'Run strict analysis' is on")
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
