#!/usr/bin/env python3
# modules.py — Per-pipeline-module config dataclasses for the runner-script GUI.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
One dataclass per pipeline module (CAAS, Disambiguation, Accumulation, RERconverge,
FADE, VEP, Scoring, Enrichment). Each field maps 1:1 to a conf/*.config param
(see the inline `# --flag` comment) so every module's full tuning surface is a real
GUI field rather than only the handful of "essential" knobs real runs vary — see
gui/models/precomputed.py's docstring for why the previously per-module standalone/
fallback ("_from"/"_input") fields moved out into one consolidated dataclass instead
of staying scattered one-per-tab.

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
    patterns: str = "1,2,3"  # --patterns
    perm_pool_size: str = "100000"  # --perm_pool_size
    caas_full_perms: str = "1000"  # --caas_full_perms
    caas_permulation_enrichment: bool = True  # --caas_permulation_enrichment

    # Contrast-selection tuning (conf/common.config) — used when --contrast_selection
    # runs upstream of CT (bundled unconditionally with CAAS, see run_single.sh.j2).
    top_quantile: str = "0.90"  # --top_quantile
    bottom_quantile: str = "0.10"  # --bottom_quantile
    contrast_max_iter: str = "3"  # --contrast_max_iter
    min_contrasts: str = "3"  # --min_contrasts

    # Discovery/resample fine-tuning (conf/ct.config)
    publish_intermediates: bool = False  # --publish_intermediates
    ct_discovery_batch_size: str = "25"  # --ct_discovery_batch_size
    ct_bootstrap_batch_size: str = "10"  # --ct_bootstrap_batch_size
    min_divergent_fraction: str = "0.5"  # --min_divergent_fraction
    max_bg_gaps_fraction: str = "0.0"  # --max_bg_gaps_fraction
    max_fg_gaps_fraction: str = "0.0"  # --max_fg_gaps_fraction
    max_gaps_fraction: str = "0.0"  # --max_gaps_fraction
    max_bg_miss_fraction: str = "0.0"  # --max_bg_miss_fraction
    max_fg_miss_fraction: str = "0.0"  # --max_fg_miss_fraction
    max_miss_fraction: str = "0.0"  # --max_miss_fraction
    miss_pair: bool = True  # --miss_pair
    caap_mode: bool = True  # --caap_mode
    fgsize: str = "6"  # --fgsize
    bgsize: str = "6"  # --bgsize
    perm_strategy: str = "BM"  # --perm_strategy (FGBG|BM|lambda)
    resample_use_n: bool = True  # --resample_use_n
    max_tries: str = "1000000"  # --max_tries
    perm_pheno_col: str = ""  # --perm_pheno_col
    traitvalues: str = ""  # --traitvalues
    chunk_size: str = "500"  # --chunk_size
    include_b0: bool = False  # --include_b0
    alpha_threshold: str = "0.05"  # --alpha_threshold

    # Debug-only (conf/ct.config warns: don't run these unless needed).
    export_groups: bool = False  # --export_groups
    export_perm_discovery: bool = False  # --export_perm_discovery

    # NOTE: discovery_from/resample_from/bootstrap_from moved to
    # gui/models/precomputed.py::PrecomputedConfig.


# ── Disambiguation (+ bundled Post-processing sub-section) ────────────────────


@dataclass(kw_only=True)
class DisambiguationConfig(ModuleConfigBase):
    ct_disambig_asr_mode: str = "precomputed"  # --ct_disambig_asr_mode (precomputed|compute)
    ct_disambig_asr_model: str = "lg"  # --ct_disambig_asr_model
    ct_disambig_asr_cache_dir: str = ""  # --ct_disambig_asr_cache_dir
    ct_disambig_convergence_mode: str = "focal_clade"  # --ct_disambig_convergence_mode (focal_clade|mrca)
    ct_disambig_posterior_threshold: str = "0.1"  # --ct_disambig_posterior_threshold
    ct_disambig_max_tasks_per_child: str = "50"  # --ct_disambig_max_tasks_per_child
    # Separate ASR Robustness diagnostics report/stage (conf/ct_disambiguation.config).
    asr_robustness: bool = True  # --asr_robustness

    # Post-processing (--ct_postproc) always runs whenever Disambiguation itself is
    # enabled — no separate toggle; see DisambiguationConfig.enabled (ModuleConfigBase).
    # Post-processing sub-section (conf/ct_postproc.config), only meaningful when enabled.
    run_postproc_exploratory: bool = True  # Run exploratory parameter sweep
    run_postproc_filter: bool = True  # Run filtering production run
    caas_postproc_mode: str = "filter"  # Legacy compatibility field
    # Filter-mode (single run) thresholds.
    filter_minlen: str = "3"  # --filter_minlen
    filter_maxcaas: str = "0.7"  # --filter_maxcaas
    # Exploratory-mode (parameter sweep) threshold lists.
    minlen_values: str = "2,3,4,10"  # --minlen_values
    maxcaas_values: str = "0.6,0.7,0.8"  # --maxcaas_values
    gene_filter_mode: str = "dubious"  # --gene_filter_mode (none|extreme|dubious|both)
    extreme_threshold: str = "0.99"  # --extreme_threshold
    iqr_multiplier: str = "3.0"  # --iqr_multiplier

    # NOTE: signification_from/disambiguation_input/disambiguation_dir/background_input
    # moved to gui/models/precomputed.py::PrecomputedConfig.


# ── Accumulation ────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class AccumulationConfig(ModuleConfigBase):
    accumulation_n_randomizations: str = "1000000"  # --accumulation_n_randomizations
    accumulation_randomization_type: str = "cons_decile"  # --accumulation_randomization_type (naive|cons_decile)
    accumulation_seed: str = "1998"  # --accumulation_seed
    accumulation_entropy_dir: str = ""  # --accumulation_entropy_dir
    accumulation_fdr: str = "0.1"  # --accumulation_fdr
    accumulation_pval_threshold: str = "0.05"  # --accumulation_pval_threshold
    accumulation_report_pval_threshold: str = "0.05"  # --accumulation_report_pval_threshold

    # NOTE: accumulation_caas_input/accumulation_background_input moved to
    # gui/models/precomputed.py::PrecomputedConfig.


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

    # Intermediate-output paths (conf/rerconverge.config) — leave blank to use the
    # pipeline's own ${params.traitname}-derived defaults.
    trait_out: str = ""  # --trait_out
    trees_out: str = ""  # --trees_out
    matrix_out: str = ""  # --matrix_out

    # Continuous-analysis tuning
    rer_minsp: str = "15"  # --rer_minsp
    winsorize_rer: str = "3"  # --winsorizeRER
    winsorize_trait: str = "3"  # --winsorizeTrait

    # Trait-type routing
    rer_trait_mode: str = "auto"  # --rer_trait_mode (auto|continuous|binary)

    # Binary-specific options
    rer_binary_clade: str = "all"  # --rer_binary_clade (all|ancestral|terminal)
    rer_min_pos: str = "2"  # --rer_min_pos

    # Report thresholds
    rer_pval_threshold: str = "0.05"  # --rer_pval_threshold
    rer_pval_column: str = "p.perm"  # --rer_pval_column
    rer_top_n_labels: str = "20"  # --rer_top_n_labels
    rer_transform: str = "ha_logit"  # --rer_transform (auto|ha_logit|logit|arcsin|log10|none)

    # RER's own tested-gene universe for the FCS report (falls back to the CAAS
    # background when unset).
    rer_universe_file: str = ""  # --rer_universe_file
    # Cross-module annotation input (SCORING's fcs_stats.tsv) for the RER FCS
    # leading-edge table.
    rer_gene_scores: str = ""  # --rer_gene_scores

    # NOTE: rer_continuous_file/rer_perms_file/scoring_rer_input/scoring_rer_perms_input
    # (used when RER itself is off) all live on gui/models/precomputed.py::PrecomputedConfig
    # (use_rer), auto-derived per phenotype from one base_path rather than typed in here
    # or per-row on PhenotypeRow.


# ── FADE ─────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class FadeConfig(ModuleConfigBase):
    # Off by default — matches RUN_FADE=false in the reference scripts.
    enabled: bool = False
    fade_mode: str = "all"  # --fade_mode (all|gene_set)

    # Shared alignment-prep / FADE-run batching (SLURM/Seqera scheduler overhead).
    selection_prep_batch_size: str = "500"  # --selection_prep_batch_size
    fade_batch_size: str = "200"  # --fade_batch_size

    # Standalone postproc gene-list paths for gene_set mode (used when FADE runs
    # without an upstream --ct_postproc step) — distinct from fade_json_dir_top/
    # bottom (a precomputed *report*, not a gene-set restriction).
    fade_postproc_top: str = ""  # --fade_postproc_top
    fade_postproc_bottom: str = ""  # --fade_postproc_bottom

    # Statistical threshold
    fade_bf_threshold: str = "100"  # --fade_bf_threshold

    # HyPhy FADE inference settings
    fade_model: str = "LG"  # --fade_model
    lg_dat_path: str = ""  # --lg_dat_path
    fade_method: str = "Variational-Bayes"  # --fade_method (Variational-Bayes|Collapsed-Gibbs|Metropolis-Hastings)
    fade_grid: str = "20"  # --fade_grid
    fade_chains: str = "5"  # --fade_chains
    fade_chain_length: str = "2000000"  # --fade_chain_length
    fade_burn_in: str = "1000000"  # --fade_burn_in
    fade_samples: str = "1000"  # --fade_samples
    fade_concentration: str = "0.5"  # --fade_concentration

    # Report options
    fade_min_genes_for_heatmap: str = "2"  # --fade_min_genes_for_heatmap
    # FADE's own tested-gene universe for the FCS report (falls back to the CAAS
    # background when unset).
    fade_universe_file: str = ""  # --fade_universe_file

    # NOTE: fade_json_dir_top/bottom/scoring_fade_summary_top/bottom/
    # scoring_fade_site_top/bottom (used when FADE itself is off) all live on
    # gui/models/precomputed.py::PrecomputedConfig (use_fade), auto-derived per
    # phenotype from one base_path rather than typed in here or per-row on
    # PhenotypeRow.


# ── VEP ──────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class VepConfig(ModuleConfigBase):
    vep_primateai_db: str = ""  # --vep_primateai_db
    vep_map_dir: str = ""  # --vep_map_dir
    # VEP's own COSMIC database (conf/vep.config) — distinct from Scoring's
    # scoring_vep_cosmic *scores* fallback (see PrecomputedConfig), which is a
    # separate param workflows/vep.nf never reads.
    cosmic_db: str = ""  # --cosmic_db

    # NOTE: vep_caas_input moved to gui/models/precomputed.py::PrecomputedConfig.


# ── Scoring ──────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ScoringConfig(ModuleConfigBase):
    scoring_stress: bool = True  # RUN_SCORING_STRESS -> --scoring_stress
    scoring_stress_top_n: str = "25"  # --scoring_stress_top_n
    scoring_stress_rank_metric: str = "spearman"  # --scoring_stress_rank_metric
    scoring_window_size_bp: str = "1000000"  # --scoring_window_size_bp
    gene_ensembl_file: str = ""  # --gene_ensembl_file
    scoring_position_top_pct: str = "0.10"  # --scoring_position_top_pct
    scoring_gene_top_pct: str = "0.10"  # --scoring_gene_top_pct

    # NOTE: scoring_postproc_input/scoring_accum_dir/scoring_vep_primateai/
    # scoring_background_input/caas_perms_file/scoring_fade_site_top/bottom moved
    # to gui/models/precomputed.py::PrecomputedConfig.


# ── Enrichment (+ bundled POSENRICH, FCS, STRING/DOMINO, COMPARE) ─────────────


@dataclass(kw_only=True)
class EnrichmentConfig(ModuleConfigBase):
    posenrich_enabled: bool = True  # RUN_POSENRICH -> --posenrich
    gmt_dir: str = ""  # --gmt_dir

    # FCS (ranked-Wilcoxon) enrichment
    fcs_min_genes: str = "5"  # --fcs_min_genes
    fcs_fdr: str = "0.15"  # --fcs_fdr
    fcs_pperm_thr: str = "0.025"  # --fcs_pperm_thr
    fcs_top_n: str = "20"  # --fcs_top_n
    # NOTE: caas_permulation_enrichment lives on CaasConfig, not here -- it's one
    # param (conf/enrichment.config) but also gates whether CT's own
    # CAAS_PERMULATION subworkflow runs (see main.nf), so its one true home is the
    # CAAS tab. A duplicate field here would silently do nothing (never wired to
    # any template flag) and imply Enrichment has independent control it doesn't.

    # STRING (ID mapping + per-DOMINO-module functional labelling only —
    # module-finding is DOMINO's job; no standalone STRING term-enrichment report)
    string_db_dir: str = ""  # --string_db_dir
    string_species: str = "9606"  # --string_species

    # DOMINO active-module identification (replaces STRING's walktrap clustering)
    domino_network_score_thr: str = "700"  # --domino_network_score_thr
    domino_slice_thr: str = "0.3"  # --domino_slice_thr
    domino_module_thr: str = "0.05"  # --domino_module_thr
    # Centralized DOMINO-based AMI run + cross-module COMPARE report
    # (conf/scoring.config, gated inside workflows/enrichment.nf). RER/FADE/
    # Accumulation's own gene lists are always computed automatically whenever
    # those tools run (no separate --ami flag) and feed this unified report's
    # FADE/RER sections directly.
    scoring_ami: bool = True  # --scoring_ami
    scoring_string: bool = True  # --scoring_string
    scoring_compare_fdr: str = "0.15"  # --scoring_compare_fdr
    scoring_compare_top_n: str = "20"  # --scoring_compare_top_n

    # POSENRICH data files (conf/enrichment.config).
    egg_members_file: str = ""  # --egg_members_file
    egg_annotations_file: str = ""  # --egg_annotations_file
    cosmic_db: str = ""  # --cosmic_db (also read by VEP; see VepConfig.cosmic_db)
    domain_variability_file: str = ""  # --domain_variability_file
    ucr_positions_file: str = ""  # --ucr_positions_file
    fubar_sites_file: str = ""  # --fubar_sites_file

    # POSENRICH thresholds (position-wise Path Sum Permulation, not the gene FCS above)
    posenrich_min_size: str = "5"  # --posenrich_min_size
    posenrich_max_size: str = "0"  # --posenrich_max_size
    posenrich_padj_thr: str = "0.15"  # --posenrich_padj_thr

    # NOTE: posenrich_background_file (CT's own background.output — used when CT
    # is off) is derived per phenotype via PrecomputedConfig.use_discovery, not a
    # separate field here. accumulation_enrichment_gene_lists_input was removed
    # entirely: no producer or consumer of it existed anywhere in the pipeline.


# ── Aggregate ─────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ModulesConfig:
    """Container for all 8 module configs, aggregated under ProjectConfig."""

    caas: CaasConfig = field(default_factory=CaasConfig)
    disambiguation: DisambiguationConfig = field(default_factory=DisambiguationConfig)
    accumulation: AccumulationConfig = field(default_factory=AccumulationConfig)
    rer: RerConfig = field(default_factory=RerConfig)
    fade: FadeConfig = field(default_factory=FadeConfig)
    vep: VepConfig = field(default_factory=VepConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
