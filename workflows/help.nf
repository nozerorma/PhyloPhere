#!/usr/bin/env nextflow

/*
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/PhyloPhere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: help.nf
#
# Usage:          nextflow run main.nf --help
#                 nextflow run main.nf --help --ct_tool discovery
#                 nextflow run main.nf --help --scoring
*/

// General help message
def general_help = '''
PHYLOPHERE Nextflow pipeline
=============================================
A complete set of phylogenetic comparative tools for Phenome-Genome studies:
convergent amino acid substitution (CAAS) discovery, ancestral-state
disambiguation, relative evolutionary rate (RERconverge), directional
selection (HyPhy FADE), VEP-style annotation, composite scoring, and
functional/network enrichment.

General usage:            > nextflow run main.nf [options]
Module-specific help:     > nextflow run main.nf --help --<module_flag>

Module flag              Description
------------------------  -------------------------------------------------
--reporting               Preliminary dataset/phenotype exploration reports
                          (independent of CAAS discovery).
--contrast_selection      Prune data and select foreground/background
                          contrasts from a continuous or discrete trait.
--ct_tool <tools>         CAAStools discovery, resample and/or bootstrap
                          (comma-separated, e.g. "discovery,resample,bootstrap").
                          See: --help --ct_tool <discovery|resample|bootstrap>
--ct_disambiguation       Classify CAAS as convergent/parallel/divergent via ASR.
--ct_postproc             Cluster/gene-level filtering + characterization report.
--ct_accumulation         Permutation test for CAAS accumulation per gene.
--vep                     Annotate CAAS positions with PrimateAI-3D / COSMIC.
--fade                    HyPhy FADE directional (top/bottom) selection scan.
--rer_tool <steps>        RERconverge relative-evolutionary-rate analysis
                          (comma-separated: build_trait,build_trees,build_matrix,continuous).
--scoring                 Composite position- and gene-level scoring across
                          CAAS, FADE, RER, accumulation and VEP signals.
--enrichment              FCS ranked enrichment, DOMINO/STRING network
                          analysis and position-level enrichment (POSENRICH).

Each module can also run standalone from a precomputed input of an upstream
module — see the module-specific help (or the README's Configuration
Reference) for the relevant "*_from" / "*_input" parameters.
'''

// ── CT: discovery / resample / bootstrap ────────────────────────────────────

def discovery_help = '''
CT Discovery — Help
=============================================
Detects Candidate Amino Acid Substitutions (CAAS) from Multiple Sequence
Alignments (MSA) against a foreground/background trait split.

Usage:
--alignment                <"input_dir">                null
--caas_config               <"caas_config_file">         null
--ali_format                <"fasta">                    "fasta"
--patterns                  <"1,2,3,4">                  "1,2,3"
--min_divergent_fraction    <FLOAT 0-1>                  0.5
--max_bg_gaps_fraction      <FLOAT 0-1>                  0.0
--max_fg_gaps_fraction      <FLOAT 0-1>                  0.0
--max_gaps_fraction         <FLOAT 0-1>                  0.0
--max_bg_miss_fraction      <FLOAT 0-1>                  0.0
--max_fg_miss_fraction      <FLOAT 0-1>                  0.0
--max_miss_fraction         <FLOAT 0-1>                  0.0
--miss_pair                 <true|false>                 true
--caap_mode                 <true|false>                 true
--ct_discovery_batch_size   <INTEGER>                    25    (genes per task)
--publish_intermediates     <true|false>                 false

NOTE: fractions are resolved to floor(n_pairs * fraction) at runtime.
'''

def resample_help = '''
CT Resample — Help
=============================================
Resamples virtual phenotypes for CAAS bootstrap analysis.

Usage:
--tree                      <"nwtree_file">              null
--perm_strategy              <"FGBG|BM|lambda">            "BM"
--resample_use_n             <true|false>                 true   (use N/C sample size count traits for Jeffreys CI filtering if available)
--fgsize                     <INTEGER>                    6
--bgsize                     <INTEGER>                    6
--traitvalues                 <"traitvalues_file">          null   (required for BM/lambda)
--perm_pheno_col             <STRING>                     ""     (trait column in --traitvalues; empty = auto-detect)
--chunk_size                  <INTEGER>                    500
--include_b0                  <true|false>                 false

Permulation sizing (three distinct quantities):
--max_tries                  <INTEGER>                    1000000 (draw budget; raised 50% up to twice if the pool is short, then fails)
--perm_pool_size             <INTEGER>                    100000  (accepted permulations harvested; the position-level null)
--caas_full_perms            <INTEGER>                    1000    (drawn from the pool for the CAAS FCS null)

Output: directory of resample_*.tab files (one per chunk_size cycles).

Strategy requirements:
FGBG                        --fgsize --bgsize
BM                           --traitvalues
'''

def bootstrap_help = '''
CT Bootstrap — Help
=============================================
Runs CAAS bootstrap analysis over resampled foreground/background labelings.

Usage:
--resample_from              <"resampleDir|resampleFile">  null
--discovery_from              <"discovery_output_dir">      null   (enables position-filter speedup)
--ct_bootstrap_batch_size     <INTEGER>                     10    (genes per task)
--export_groups               <true|false>                  false  (debug)
--export_perm_discovery       <true|false>                  false  (debug; required for CAAS permulation-excess null)
--alpha_threshold             <FLOAT>                        0.05

# Shares the same alignment/filter parameters as discovery:
--alignment, --caas_config, --ali_format, --patterns,
--min_divergent_fraction, --max_bg_gaps_fraction, --max_fg_gaps_fraction,
--max_gaps_fraction, --max_bg_miss_fraction, --max_fg_miss_fraction,
--max_miss_fraction, --miss_pair, --caap_mode
'''

// ── Contrast selection ──────────────────────────────────────────────────────

def contrast_selection_help = '''
Contrast Selection — Help
=============================================
Prunes the trait/tree pair and selects foreground/background extremes via
an independent-contrasts algorithm before CAAS discovery.

Usage:
--my_traits                  <"traits.csv">               null
--traitname                   <"column_name">               null
--sp_colname                  <"species">                    "species"
--trait_type                   <"auto|continuous|ordinal">   "auto"  (auto-infers; continuous = PSS, ordinal = coded fg/bg)
--pss_top_pct                   <FLOAT 0-1>                   0.05   (top PSS percentile gate, continuous traits)
--perm_strategy                  <"auto|OU|BM">               "auto"  (evolutionary model for PSS / permulation)
--max_contrasts                   <INTEGER>                     0     (0 = dynamic discovery)
--min_contrasts                    <INTEGER>                     3     (min foreground pairs required)
--prune_data                        <true|false>                   false  (requires --reporting)
--prune_list                         <"species_list">               null
--prune_list_secondary                 <"species_list">               null
'''

// ── CT disambiguation / post-processing / accumulation ─────────────────────

def disambiguation_help = '''
CT Disambiguation — Help
=============================================
Classifies each discovered CAAS as convergent, parallel or divergent using
ancestral state reconstruction (ASR) relative to the phenotype tree topology.

Usage:
--signification_from            <"signification_output">     null  (standalone entry)
--ct_disambig_asr_mode            <"precomputed">                "precomputed"
--ct_disambig_asr_model             <"lg">                          "lg"
--ct_disambig_asr_cache_dir          <"cache_dir">                   null
--ct_disambig_convergence_mode        <"focal_clade|mrca">            "focal_clade"
--ct_disambig_posterior_threshold      <FLOAT 0-1>                     0.1
--ct_disambig_max_tasks_per_child       <INTEGER>                       50
--asr_robustness                        <true|false>                    true   (parallel diagnostic report)
'''

def postproc_help = '''
CT Post-processing — Help
=============================================
Cluster/gene-level filtering of disambiguated CAAS, background cleanup, and
the Characterization report.

Usage:
--disambiguation_input        <"caas_convergence_master.csv">   null  (standalone entry)
--disambiguation_dir            <"disambiguation_dir|.tar.gz">     null
--background_input                <"background_file">                 null
--caas_postproc_mode                <"exploratory|filter">              "filter"
--minlen_values                       <"2,3,4,10">                        (exploratory mode)
--maxcaas_values                        <"0.6,0.7,0.8">                     (exploratory mode)
--filter_minlen                          <INTEGER>                            3     (filter mode)
--filter_maxcaas                          <FLOAT 0-1>                          0.7   (filter mode)
--gene_filter_mode                          <"none|extreme|dubious|both">        "dubious"
--extreme_threshold                           <FLOAT 0-1>                          0.99
--iqr_multiplier                                <FLOAT>                              3.0
'''

def accumulation_help = '''
CT Accumulation — Help
=============================================
Permutation test for whether a gene accumulates more CAAS than expected by
chance; runs once per direction (top/bottom/all).

Usage:
--accumulation_caas_input             <"filtered_discovery.tsv">   null  (standalone entry)
--accumulation_background_input        <"background_file">           null
--accumulation_entropy_dir               <"entropy_dir">                null
--accumulation_randomization_type          <"naive|cons_decile">          "cons_decile"
--accumulation_n_randomizations              <INTEGER>                      1000000
--accumulation_fdr                               <FLOAT 0-1>                    0.1
'''

// ── VEP ──────────────────────────────────────────────────────────────────────

def vep_help = '''
VEP Characterization — Help
=============================================
Annotates filtered CAAS positions with PrimateAI-3D pathogenicity scores and
COSMIC somatic-mutation evidence (each skipped if its DB param is empty).

Usage:
--vep_caas_input               <"caas_input">          null  (standalone entry)
--vep_map_dir                    <"map_dir">              null  (required when --vep is set)
--vep_primateai_db                 <"primateai_db">          null
--cosmic_db                          <"cosmic_db">              null
'''

// ── FADE ─────────────────────────────────────────────────────────────────────

def fade_help = '''
FADE (directional selection) — Help
=============================================
Runs HyPhy FADE, a Bayesian branch-site model, to detect accelerated or
decelerated amino-acid selection on phenotype-extreme branches. Always runs
both "top" and "bottom" directions.

Usage:
--selection_prep_batch_size          <INTEGER>                500
--fade_batch_size                      <INTEGER>                200
--fade_bf_threshold                      <INTEGER>                100  (Bayes Factor cutoff)
--fade_model                               <"LG">                    "LG"
--fade_method                                <"Variational-Bayes|Collapsed-Gibbs|Metropolis-Hastings">  "Variational-Bayes"
--fade_grid                                    <INTEGER>                20
--fade_chains, --fade_chain_length,
--fade_burn_in, --fade_samples               (MCMC-only settings; ignored by Variational-Bayes)
--fade_concentration                            <FLOAT>                  0.5
--fade_universe_file                              <"gene_list">            null  (else falls back to cleaned background)
--fade_json_dir_top / --fade_json_dir_bottom        <"dir/of/*.FADE.json">   null  (render report from precomputed JSON)
'''

// ── RERconverge ──────────────────────────────────────────────────────────────

def rer_help = '''
RERconverge (RER) — Help
=============================================
Computes Relative Evolutionary Rate — branch-length deviation correlated with
the trait — auto-routing to continuous or binary correlation.

Usage:
--rer_tool                    <"build_trait,build_trees,build_matrix,continuous">   ""
--rer_trait_mode                 <"auto|continuous|binary">                          "auto"
--rer_minsp                        <INTEGER>                                            15
--winsorizeRER, --winsorizeTrait     <INTEGER>                                            3
--rer_perm_batches                    <INTEGER>                                            10   (0 = skip permulation)
--rer_perms_per_batch                   <INTEGER>                                            100
--rer_binary_clade                        <"all|ancestral|terminal">                           "all"
--rer_min_pos                               <INTEGER>                                            2
--rer_pval_threshold                          <FLOAT 0-1>                                          0.05
--rer_transform                                 <"auto|ha_logit|logit|arcsin|log10|none">            "ha_logit"
--rer_universe_file                               <"gene_list">                                        null  (else falls back to cleaned background)
--rer_continuous_file                               <"trait.continuous.output">                          null  (render report from precomputed RDS)
--rer_perms_file                                      <"trait.continuous.perms.rds">                        null
'''

// ── Scoring ──────────────────────────────────────────────────────────────────

def scoring_help = """
CAAS Scoring — Help
=============================================
Computes composite CAAS scores at position-level and gene-level, integrating
CT, FADE, RER, accumulation and VEP signals.

Usage:
--scoring                   <true|false>            false
--scoring_ami                <true|false>            true    (downstream DOMINO/AMI active module analysis)
--scoring_string             <true|false>            true    (legacy alias for --scoring_ami)
--scoring_compare_fdr          <FLOAT 0-1>             0.15
--scoring_compare_top_n           <INTEGER>               20
--scoring_stress                   <true|false>            true
--scoring_stress_top_n                <INTEGER>               25
--scoring_stress_rank_metric            <"spearman|pearson">    "spearman"
--scoring_position_top_pct                <FLOAT 0-1>             0.10
--scoring_gene_top_pct                       <FLOAT 0-1>             0.10
--scoring_gene_perm_pooled                <true|false>            false  (opt-in n-stratified pooled-null variant of gene_caas_pperm, higher resolution than the per-gene-row floor -- valid only within an n_positions stratum)
--scoring_window_size_bp                       <INTEGER>               1000000

Standalone mode (provide inputs directly):
--scoring_postproc_input        <"filtered_discovery.tsv">        ""
--scoring_background_input        <"background.txt">                ""
--scoring_fade_summary_top          <"fade_summary_top.tsv">           ""
--scoring_fade_summary_bottom          <"fade_summary_bottom.tsv">        ""
--scoring_fade_site_top                  <"fade_site_top.tsv">              ""
--scoring_fade_site_bottom                  <"fade_site_bottom.tsv">           ""
--scoring_rer_input                           <"rerconverge_summary.tsv">        ""
--scoring_accum_dir                              <"accumulation_dir/">              ""
--scoring_vep_primateai                            <"vep_primateai.tsv">              ""    (also POSENRICH's precomputed-VEP fallback, see enrichment.nf)
--scoring_vep_cosmic                                  <"vep_cosmic.tsv">                  ""    (also POSENRICH's precomputed-VEP fallback, see enrichment.nf)
--caas_perms_file                                        <"caas_perms.rds">                  ""    (genome-wide CAAS permulation null)
--caas_pos_detail_file                                   <"perm_pos_detail/">          ""    (rebuild the null from a prior run's per-cycle detail shard dir, or a legacy perm_pos_detail.tsv.gz; preferred over --caas_perms_file for reruns, minutes, no ASR replay)

Position-level components: biochem, ASR, convergence, parallel, [FADE]
Gene-level scores: gene_caas, [gene_rer], [gene_fade], gene_composite
"""

// ── Enrichment ───────────────────────────────────────────────────────────────

def enrichment_help = '''
Enrichment (FCS / DOMINO-STRING / POSENRICH) — Help
=============================================
Ranked gene-set enrichment (FCS, Wilcoxon-AUC on GMT gene sets) run
separately for CAAS and RER rankings, an orthogonal cross-module Comparison
report, optional DOMINO active-module detection + STRING functional
enrichment, and optional position-level Fisher-exact enrichment (POSENRICH).

Usage:
--enrichment                        <true|false>          false
--gmt_dir                            <"gmt_dir">              null
--fcs_min_genes                        <INTEGER>                5
--fcs_max_genes                        <INTEGER>                500   (0 = no cap)
--fcs_fdr                                <FLOAT 0-1>              0.15
--fcs_pperm_thr                            <FLOAT 0-1>              0.025
--fcs_top_n                                  <INTEGER>                20
--caas_permulation_enrichment                    <true|false>             true

DOMINO / STRING:
--string_species                    <NCBI taxon id>          9606
--string_db_dir                           <"string_cache_dir">     null
--domino_network_score_thr                  <INTEGER>                700
--domino_slice_thr                             <FLOAT 0-1>              0.3
--domino_module_thr                              <FLOAT 0-1>              0.05

POSENRICH (position-level):
--posenrich                        <true|false>          true
--posenrich_min_size                 <INTEGER>                5
--posenrich_max_size                   <INTEGER>                0    (0 = no cap)
--posenrich_padj_thr                     <FLOAT 0-1>              0.15
--posenrich_background_file                  <"background_file">      null
--domain_variability_file, --ucr_positions_file, --fubar_sites_file,
--egg_members_file, --egg_annotations_file   (position-level annotation sources)
'''

workflow HELP {
    // Check if --help is provided
    if (params.help) {
        // Most module toggles default to `true`/non-empty in conf/common.config and
        // conf/ct.config (the "core" CT+disambiguation+postproc+accumulation+VEP chain
        // runs unless explicitly disabled), so there is no boolean signal that reliably
        // means "the user explicitly asked for this module." Instead: always show the
        // general overview first, then append detail for every module that WOULD run
        // given the current params — i.e. help for the run you are about to launch.
        // To see only the overview, explicitly disable the modules you don't need
        // (e.g. --ct_tool "" --ct_disambiguation false --ct_postproc false
        // --ct_accumulation false --vep false --help).
        log.info general_help

        if (params.contrast_selection)   log.info contrast_selection_help

        // --ct_tool accepts a comma-separated list (e.g. "discovery,resample,bootstrap")
        if (params.ct_tool) {
            def tools = params.ct_tool.toString().split(',').collect { it.trim() }
            if ('discovery' in tools) log.info discovery_help
            if ('resample'  in tools) log.info resample_help
            if ('bootstrap' in tools) log.info bootstrap_help
        }

        if (params.ct_disambiguation)    log.info disambiguation_help
        if (params.ct_postproc)          log.info postproc_help
        if (params.ct_accumulation)      log.info accumulation_help
        if (params.vep)                  log.info vep_help
        if (params.fade)                 log.info fade_help
        if (params.rer_tool)             log.info rer_help
        if (params.scoring)              log.info scoring_help
        if (params.enrichment)           log.info enrichment_help

        exit 1
    }
}
