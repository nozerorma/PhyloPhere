#!/usr/bin/env python3
# modules.py — Per-pipeline-module config dataclasses for the runner-script GUI.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
One dataclass per pipeline module (CAAS, Disambiguation, Accumulation, RERconverge,
FADE, VEP, Scoring, Enrichment). Each holds only the "essential" fields actually
varied across real runs (see SBATCH_run_phenotypes_primates.sh / run_phenotype_single_primates.sh)
plus the standalone `_from`/`_input` fallback paths a downstream module needs when
this module is disabled (see conf/*.config), and a free-text `extra_flags` escape
hatch for anything not exposed as a dedicated field.

No PySide6 imports here — this module must stay importable headless (see
gui/generation/, which renders these into shell scripts without any GUI dependency).
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field


# ── Base ──────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ModuleConfigBase:
    """Shared shape for every module tab: an enable toggle and a raw-flags override."""

    enabled: bool = True
    # Free-text power-user override, appended verbatim to this module's NF_FLAGS block.
    extra_flags: str = ""


# ── CAAS / CT (contrast selection + discovery/resample/bootstrap) ─────────────


@dataclass(kw_only=True)
class CaasConfig(ModuleConfigBase):
    ct_tool_discovery: bool = True
    ct_tool_resample: bool = True
    ct_tool_bootstrap: bool = True
    caas_config_path: str = ""  # --caas_config
    perms_cycles: str = "1000000"  # --perms_cycles
    caas_full_perms: str = "1000"  # --caas_full_perms
    caas_permulation_enrichment: bool = True  # --caas_permulation_enrichment

    # Fallback inputs for downstream modules when this one is disabled.
    discovery_from: str = ""  # --discovery_from
    resample_from: str = ""  # --resample_from
    bootstrap_from: str = ""  # --bootstrap_from


# ── Disambiguation (+ bundled Post-processing sub-section) ────────────────────


@dataclass(kw_only=True)
class DisambiguationConfig(ModuleConfigBase):
    ct_disambig_asr_mode: str = "precomputed"  # --ct_disambig_asr_mode
    ct_disambig_asr_cache_dir: str = ""  # --ct_disambig_asr_cache_dir
    ct_postproc: bool = True  # --ct_postproc (bundled toggle, matches reference script)

    # Post-processing sub-section (conf/ct_postproc.config), only meaningful if ct_postproc=True.
    caas_postproc_mode: str = "filter"  # --caas_postproc_mode (filter|exploratory)
    filter_minlen: str = "3"  # --filter_minlen
    filter_maxcaas: str = "0.7"  # --filter_maxcaas
    gene_filter_mode: str = "dubious"  # --gene_filter_mode

    # Fallback inputs for downstream modules when this one is disabled.
    signification_from: str = ""  # --signification_from
    disambiguation_input: str = ""  # --disambiguation_input
    disambiguation_dir: str = ""  # --disambiguation_dir
    background_input: str = ""  # --background_input


# ── Accumulation ────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class AccumulationConfig(ModuleConfigBase):
    accumulation_n_randomizations: str = "1000000"  # --accumulation_n_randomizations
    accumulation_entropy_dir: str = ""  # --accumulation_entropy_dir

    # Fallback inputs for downstream modules when this one is disabled.
    accumulation_caas_input: str = ""  # --accumulation_caas_input
    accumulation_background_input: str = ""  # --accumulation_background_input


# ── RERconverge ─────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class RerConfig(ModuleConfigBase):
    # Off by default — matches RUN_RER=false in the reference scripts.
    enabled: bool = False
    rer_tool_build_trait: bool = True
    rer_tool_build_tree: bool = True
    rer_tool_build_matrix: bool = True
    rer_tool_continuous: bool = True
    gene_trees: str = ""  # --gene_trees
    rer_perm_batches: str = "10"  # --rer_perm_batches
    rer_perms_per_batch: str = "100"  # --rer_perms_per_batch

    # NOTE: no scoring_rer_input here — the reference script sets it per phenotype
    # (see gui/models/runtime.py's PhenotypeRow), not once globally for the module.


# ── FADE ─────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class FadeConfig(ModuleConfigBase):
    # Off by default — matches RUN_FADE=false in the reference scripts.
    enabled: bool = False
    fade_mode: str = "all"  # --fade_mode (all|gene_set)

    # NOTE: no scoring_fade_summary_top/bottom here — the reference script sets
    # those per phenotype (see gui/models/runtime.py's PhenotypeRow). The
    # non-per-row scoring_fade_site_top/bottom fallback lives on ScoringConfig.


# ── VEP ──────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class VepConfig(ModuleConfigBase):
    vep_primateai_db: str = ""  # --vep_primateai_db
    vep_map_dir: str = ""  # --vep_map_dir
    scoring_vep_cosmic: str = ""  # --scoring_vep_cosmic (COSMIC db, threaded through scoring)

    # Fallback input for Scoring when VEP is disabled.
    vep_caas_input: str = ""  # --vep_caas_input


# ── Scoring ──────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ScoringConfig(ModuleConfigBase):
    scoring_stress: bool = True  # RUN_SCORING_STRESS -> --scoring_stress
    scoring_window_size_bp: str = "1000000"  # --scoring_window_size_bp
    gene_ensembl_file: str = ""  # --gene_ensembl_file

    # Fallback inputs when Scoring runs standalone against precomputed upstream outputs.
    scoring_postproc_input: str = ""  # --scoring_postproc_input
    scoring_accum_dir: str = ""  # --scoring_accum_dir
    scoring_vep_primateai: str = ""  # --scoring_vep_primateai
    scoring_background_input: str = ""  # --scoring_background_input
    caas_perms_file: str = ""  # --caas_perms_file
    # Global (non-per-phenotype) FADE site-level fallback, unlike scoring_fade_summary_top/
    # bottom which are set per phenotype row (see PhenotypeRow).
    scoring_fade_site_top: str = ""  # --scoring_fade_site_top
    scoring_fade_site_bottom: str = ""  # --scoring_fade_site_bottom


# ── Enrichment (+ bundled POSENRICH) ──────────────────────────────────────────


@dataclass(kw_only=True)
class EnrichmentConfig(ModuleConfigBase):
    posenrich_enabled: bool = True  # RUN_POSENRICH -> --posenrich
    gmt_dir: str = ""  # --gmt_dir
    string_db_dir: str = ""  # --string_db_dir

    # POSENRICH data files (conf/enrichment.config).
    egg_members_file: str = ""  # --egg_members_file
    egg_annotations_file: str = ""  # --egg_annotations_file
    cosmic_db: str = ""  # --cosmic_db
    domain_variability_file: str = ""  # --domain_variability_file
    ucr_positions_file: str = ""  # --ucr_positions_file
    fubar_sites_file: str = ""  # --fubar_sites_file

    # Fallback inputs when Enrichment runs standalone.
    enrichment_background_input: str = ""  # --enrichment_background_input
    enrichment_gene_lists_input: str = ""  # --enrichment_gene_lists_input
    accumulation_enrichment_gene_lists_input: str = ""  # --accumulation_enrichment_gene_lists_input
    posenrich_background_file: str = ""  # --posenrich_background_file


# ── Aggregate ─────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ModulesConfig:
    """Container for all 9 module configs, aggregated under ProjectConfig."""

    caas: CaasConfig = field(default_factory=CaasConfig)
    disambiguation: DisambiguationConfig = field(default_factory=DisambiguationConfig)
    accumulation: AccumulationConfig = field(default_factory=AccumulationConfig)
    rer: RerConfig = field(default_factory=RerConfig)
    fade: FadeConfig = field(default_factory=FadeConfig)
    vep: VepConfig = field(default_factory=VepConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
