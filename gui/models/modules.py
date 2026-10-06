#!/usr/bin/env python3
# modules.py — Per-pipeline-module config dataclasses for the runner-script GUI.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Modules: one dataclass per pipeline module (CAAS, Disambiguation, Accumulation,
RERconverge, FADE, VEP, Scoring, Enrichment) plus the ModulesConfig container.

Each field maps to one conf/*.config param (named in the inline `# --flag` comment),
so the full tuning surface of every module is a GUI field. Inputs that substitute for
an upstream module's output (the "_from" and "_input" params) are not fields here:
they are derived per phenotype from PrecomputedConfig (gui/models/precomputed.py).

No PySide6 import: the module stays importable headless, because gui/generation/
renders these dataclasses into shell scripts without a GUI.

Imported by: gui/models/project.py, gui/widgets/common/module_tab.py, the module tabs
in gui/widgets/tabs/
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field


# ── Base ──────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ModuleConfigBase:
    """Shared shape of every module tab: an enable toggle and a raw-flags override."""

    enabled: bool = True
    # Free text appended verbatim to the nextflow command line of the generated script
    # (run_single.sh.j2), when this module is enabled.
    extra_flags: str = ""


# ── CAAS / CT (contrast selection + discovery/resample) ───────────────────────


@dataclass(kw_only=True)
class CaasConfig(ModuleConfigBase):
    ct_tool_discovery: bool = True
    ct_tool_resample: bool = True
    caas_config_path: str = ""  # --caas_config
    patterns: str = "1,2,3"  # --patterns
    caas_full_perms: str = "1000"  # --caas_full_perms
    caas_permulation_enrichment: bool = True  # --caas_permulation_enrichment
    caas_perms_postproc: bool = True  # --caas_perms_postproc (apply the observed post-processing filters to the permulation null)

    # Contrast-selection tuning (conf/common.config), used when --contrast_selection
    # runs upstream of CT.
    pss_top_pct: str = "0.05"  # --pss_top_pct
    max_contrasts: str = "0"  # --max_contrasts (0 = dynamic discovery)
    min_contrasts: str = "3"  # --min_contrasts

    # Discovery/resample fine-tuning (conf/ct.config)
    publish_intermediates: bool = False  # --publish_intermediates
    ct_core_batch_size: str = "20"  # --ct_core_batch_size (genes per CAAS_CORE_BATCHED task of the permulation null; 1 = one task per gene)
    min_divergent_fraction: str = "0.5"  # --min_divergent_fraction
    max_bg_gaps_fraction: str = "0.0"  # --max_bg_gaps_fraction
    max_fg_gaps_fraction: str = "0.0"  # --max_fg_gaps_fraction
    max_gaps_fraction: str = "0.0"  # --max_gaps_fraction
    max_bg_miss_fraction: str = "0.0"  # --max_bg_miss_fraction
    max_fg_miss_fraction: str = "0.0"  # --max_fg_miss_fraction
    max_miss_fraction: str = "0.0"  # --max_miss_fraction
    miss_pair: bool = True  # --miss_pair
    caap_mode: bool = True  # --caap_mode
    perm_strategy: str = "auto"  # --perm_strategy (auto|OU|BM)
    max_tries: str = "1000000"  # --max_tries
    chunk_size: str = "500"  # --chunk_size
    resample_use_n: bool = True  # --resample_use_n
    multi_hypothesis: bool = True  # --multi_hypothesis
    max_fop: str = "100"  # --max_fop (max FOP alternative hypotheses H1..Hn per contrast)

    # Precomputed discovery/resample inputs (discovery_from, resample_from) are
    # derived from PrecomputedConfig.


# ── Disambiguation (+ bundled Post-processing sub-section) ────────────────────


@dataclass(kw_only=True)
class DisambiguationConfig(ModuleConfigBase):
    ct_disambig_asr_model: str = "lg"  # --ct_disambig_asr_model
    ct_disambig_asr_cache_dir: str = ""  # --ct_disambig_asr_cache_dir
    ct_disambig_hypotheses_pairs: str = ""  # --ct_disambig_hypotheses_pairs (contrast_hypotheses_pairs.tsv override for the observed scoring)
    ct_disambig_posterior_threshold: str = "0.1"  # --ct_disambig_posterior_threshold
    ct_disambig_max_tasks_per_child: str = "50"  # --ct_disambig_max_tasks_per_child
    # ASR robustness diagnostics report, a separate stage (conf/ct_disambiguation.config).
    asr_robustness: bool = True  # --asr_robustness

    # Post-processing (--ct_postproc) has no toggle of its own: it runs whenever
    # Disambiguation is enabled, unless PrecomputedConfig.use_postproc supplies its
    # outputs (see gui/generation/context.py). Parameters of conf/ct_postproc.config.
    run_postproc_exploratory: bool = True  # Run exploratory parameter sweep
    run_postproc_filter: bool = True  # Run filtering production run
    caas_postproc_mode: str = "filter"  # --caas_postproc_mode (filter|exploratory); the generated scripts set it per pass
    # Filter-mode (single run) thresholds.
    filter_minlen: str = "3"  # --filter_minlen
    filter_maxcaas: str = "0.7"  # --filter_maxcaas
    # Directory of the trimmer's per-gene MAP tables (--caas_map_dir). Post-processing
    # uses it to place clusters in untrimmed alignment columns, for the observed and the null
    # sets; VEP reads it from here as well.
    caas_map_dir: str = ""  # --caas_map_dir
    # Exploratory-mode (parameter sweep) threshold lists.
    minlen_values: str = "2,3,4,10"  # --minlen_values
    maxcaas_values: str = "0.6,0.7,0.8"  # --maxcaas_values
    gene_filter_mode: str = "dubious"  # --gene_filter_mode (none|extreme|dubious|both)
    remove_caas_clusters: bool = True  # --remove_caas_clusters (discard spatial clusters)
    extreme_threshold: str = "0.99"  # --extreme_threshold
    iqr_multiplier: str = "3.0"  # --iqr_multiplier

    # Precomputed inputs (meta_caas_from, disambiguation_input, disambiguation_dir,
    # background_input) are derived from PrecomputedConfig.


# ── Accumulation ──────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class AccumulationConfig(ModuleConfigBase):
    accumulation_n_randomizations: str = "1000000"  # --accumulation_n_randomizations
    accumulation_randomization_type: str = "cons_decile"  # --accumulation_randomization_type (naive|cons_decile|permulation)
    accumulation_entropy_dir: str = ""  # --accumulation_entropy_dir
    accumulation_fdr: str = "0.1"  # --accumulation_fdr

    # Precomputed inputs (accumulation_caas_input, accumulation_background_input) are
    # derived from PrecomputedConfig.


# ── RERconverge ───────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class RerConfig(ModuleConfigBase):
    # Off by default, as RUN_RER=false in example_sbatch_run_phenotypes.sh.
    enabled: bool = False
    rer_tool_build_trait: bool = True
    rer_tool_build_tree: bool = True
    rer_tool_build_matrix: bool = True
    rer_tool_continuous: bool = True
    gene_trees: str = ""  # --gene_trees
    rer_perm_batches: str = "10"  # --rer_perm_batches
    rer_perms_per_batch: str = "100"  # --rer_perms_per_batch

    # Intermediate-output paths (conf/rerconverge.config). Blank selects the
    # defaults derived from ${params.traitname}.
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

    # Gene universe of the RER FCS report. Unset, the genes RERconverge tested in
    # this run are used, or the CAAS background when RER did not run (workflows/enrichment.nf).
    rer_universe_file: str = ""  # --rer_universe_file
    # SCORING's fcs_stats.tsv, used for the cross-module flags of the RER FCS
    # leading-edge table.
    rer_gene_scores: str = ""  # --rer_gene_scores

    # Precomputed RER inputs (rer_continuous_file, rer_perms_file, scoring_rer_input,
    # scoring_rer_perms_input) are derived from PrecomputedConfig (use_rer).


# ── FADE ──────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class FadeConfig(ModuleConfigBase):
    # Off by default, as RUN_FADE=false in example_sbatch_run_phenotypes.sh.
    enabled: bool = False

    # Direction and background scope
    fade_direction: str = "both"  # --fade_direction (top|bottom|both)
    fade_background_scope: str = "all"  # --fade_background_scope (all|opposite)
    fade_internal_nodes: str = "all_descendants"  # --fade_internal_nodes (all_descendants|none)
    fade_species_file: str = ""  # --fade_species_file (optional user-supplied fg/bg species file, candidate_species.tab format)

    # Batch sizes of the shared alignment-preparation and FADE-run tasks (genes per
    # task; larger batches reduce scheduler overhead).
    selection_prep_batch_size: str = "500"  # --selection_prep_batch_size
    fade_batch_size: str = "200"  # --fade_batch_size

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
    # Gene universe of the FCS report for FADE. Unset, the union of the tested genes
    # of both directions is used (workflows/enrichment.nf).
    fade_universe_file: str = ""  # --fade_universe_file

    # Precomputed FADE inputs (fade_json_dir_top/bottom, scoring_fade_summary_top/bottom,
    # scoring_fade_site_top/bottom) are derived from PrecomputedConfig (use_fade).


# ── VEP ───────────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class VepConfig(ModuleConfigBase):
    vep_primateai_db: str = ""  # --vep_primateai_db
    # The per-gene MAP directory is DisambiguationConfig.caas_map_dir (--caas_map_dir);
    # VEP reads it from there.
    # COSMIC mutation database that VEP maps positions onto (conf/vep.config). It is not
    # the scoring_vep_cosmic scores table that SCORING can take as a precomputed input
    # (see PrecomputedConfig), which workflows/vep.nf does not read.
    cosmic_db: str = ""  # --cosmic_db
    vep_ensembl: bool = False  # --vep_ensembl
    vep_cache_dir: str = ""  # --vep_cache_dir
    vep_species: str = "homo_sapiens"  # --vep_species
    vep_assembly: str = "GRCh38"  # --vep_assembly



# ── Scoring ───────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ScoringConfig(ModuleConfigBase):
    scoring_gene_top_pct: str = "0.10"  # --scoring_gene_top_pct
    scoring_position_top_pct: str = "0.10"  # --scoring_position_top_pct
    gene_ensembl_file: str = ""  # --gene_ensembl_file
    auto_generate_ensembl: bool = False  # Generate gene_ensembl_file via BioMart if unset
    ensembl_dataset: str = ""  # Ensembl BioMart dataset (e.g. hsapiens_gene_ensembl); blank -> derived from ref_species_name

    # Advanced parameters (conf/scoring.config)
    scoring_window_size_bp: str = "1000000"  # --scoring_window_size_bp
    scoring_p_emp_thr: str = "0.05"  # --scoring_p_emp_thr (position-level CAAS permulation p.adj_bh / p.adj_sam; also gates gene_caas_pperm_adj)
    caas_evidence_top_n: str = "0"  # --caas_evidence_top_n (evidence table of the N best positions after SCORING; 0 = off)
    scoring_gene_perm_pooled: bool = False  # --scoring_gene_perm_pooled (opt-in n-stratified pooled-null gene permulation p)
    scoring_hypotheses_pairs: str = ""  # --scoring_hypotheses_pairs (contrast_hypotheses_pairs.tsv override for SCORING)

    # Precomputed inputs (scoring_postproc_input, scoring_accum_dir, scoring_vep_primateai,
    # scoring_background_input, caas_perms_file, scoring_fade_site_top/bottom) are derived
    # from PrecomputedConfig.


# ── Enrichment (+ bundled POSENRICH, FCS, STRING/DOMINO, COMPARE) ─────────────


@dataclass(kw_only=True)
class EnrichmentConfig(ModuleConfigBase):
    posenrich_enabled: bool = True  # RUN_POSENRICH -> --posenrich
    fcs_enabled: bool = True  # Gate FCS (ranked-Wilcoxon) enrichment
    gmt_dir: str = ""  # --gmt_dir
    auto_fetch_gmt: bool = False  # Also download the current GO and WikiPathways gene sets
    auto_fetch_eggnog: bool = False  # Download the eggNOG pair instead of using the versioned copy

    # FCS (ranked-Wilcoxon) enrichment
    fcs_min_genes: str = "5"  # --fcs_min_genes
    fcs_max_genes: str = "1000"  # --fcs_max_genes (0 = no limit)
    fcs_fdr: str = "0.15"  # --fcs_fdr
    fcs_fdr_wilcoxon: str = ""  # --fcs_fdr_wilcoxon (blank -> the FCS FDR)
    fcs_fdr_lachenbruch: str = ""  # --fcs_fdr_lachenbruch (blank -> the FCS FDR)
    fcs_fdr_permsum: str = "0.05"  # --fcs_fdr_permsum
    pfam_cache_dir: str = ""  # --pfam_cache_dir (blank -> ~/.cache/phylophere/pfam)
    fcs_pperm_thr: str = "0.025"  # --fcs_pperm_thr
    fcs_top_n: str = "20"  # --fcs_top_n
    fcs_batch_size: str = "4"  # --fcs_batch_size (GMTs per FCS_COMPUTE_BATCHED task)
    # caas_permulation_enrichment (conf/enrichment.config) is a field of CaasConfig, not
    # of this class: it also gates whether CT's permulation core (CAAS_CORE) runs
    # (see main.nf), so the CAAS tab owns it.

    # STRING: ID mapping and functional labelling of each DOMINO module. Module
    # detection itself is done by DOMINO.
    string_db_dir: str = ""  # --string_db_dir
    string_cache_dir: str = ""  # --string_cache_dir
    string_species: str = ""  # --string_species (blank -> fallback to RuntimeConfig.ref_species_taxid)

    # DOMINO active-module identification
    domino_network_score_thr: str = "700"  # --domino_network_score_thr
    domino_slice_thr: str = "0.3"  # --domino_slice_thr
    domino_module_thr: str = "0.05"  # --domino_module_thr
    # Off by default: the network, module and edge-score files are large and only
    # needed to regenerate the AMI report by hand.
    publish_domino_intermediates: bool = False  # --publish_domino_intermediates
    # Active-module identification (DOMINO) run and cross-module COMPARE report
    # (conf/scoring.config, gated inside workflows/enrichment.nf). The gene lists of
    # RER, FADE and Accumulation are built whenever those tools run and feed the
    # FADE and RER sections of this report directly.
    scoring_ami: bool = True  # --scoring_ami
    scoring_string: bool = True  # --scoring_string
    scoring_compare_fdr: str = "0.15"  # --scoring_compare_fdr
    scoring_compare_top_n: str = "20"  # --scoring_compare_top_n
    # Concordance null of the unpaired CAAS and RER gene scores in the COMPARE
    # report (a randomization null, not a joint permulation p).
    comparison_perm_null: bool = True  # --comparison_perm_null
    comparison_perm_stat: str = "spearman"  # --comparison_perm_stat ("spearman" | "topk_overlap")
    comparison_perm_topk: str = "0.05"  # --comparison_perm_topk

    # POSENRICH component toggles & data files (conf/enrichment.config).
    posenrich_domains: bool = True  # Run Pfam domain variability analysis
    domain_variability_file: str = ""  # --domain_variability_file
    domain_ref_species: str = ""  # --domain_ref_species (blank -> fallback to RuntimeConfig.ref_species_name)

    posenrich_ucr: bool = True  # Run ultraconserved regions (UCR) analysis
    ucr_positions_file: str = ""  # --ucr_positions_file

    posenrich_eggnog: bool = False  # Run eggNOG ortholog position mapping
    eggnog_taxid: str = ""  # eggNOG taxon level (blank -> fallback to RuntimeConfig.clade_taxid)
    egg_members_file: str = ""  # --egg_members_file
    egg_annotations_file: str = ""  # --egg_annotations_file

    fubar_sites_file: str = ""  # --fubar_sites_file

    # POSENRICH thresholds (position-wise Path Sum Permulation, not the gene FCS above)
    posenrich_min_size: str = "5"  # --posenrich_min_size
    posenrich_max_size: str = "0"  # --posenrich_max_size
    posenrich_padj_thr: str = "0.05"  # --posenrich_padj_thr (Permsum-family test: no independent analytic estimate to dual-gate against, same 0.05 fdr_permsum uses in FCS)
    posenrich_batch_size: str = "4"  # --posenrich_batch_size (GMTs per POSENRICH_RUN_BATCHED task; 1 = no batching)

    # posenrich_background_file (CT's background.output, used when CT is off) is derived
    # from PrecomputedConfig.use_discovery.


# ── Aggregate ─────────────────────────────────────────────────────────────────


@dataclass(kw_only=True)
class ModulesConfig:
    """The 8 module configs, held by ProjectConfig.modules."""

    caas: CaasConfig = field(default_factory=CaasConfig)
    disambiguation: DisambiguationConfig = field(default_factory=DisambiguationConfig)
    accumulation: AccumulationConfig = field(default_factory=AccumulationConfig)
    rer: RerConfig = field(default_factory=RerConfig)
    fade: FadeConfig = field(default_factory=FadeConfig)
    vep: VepConfig = field(default_factory=VepConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)
