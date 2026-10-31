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
def general_help() {
    return '''
PHYLOPHERE Nextflow pipeline
=============================================
A complete set of phylogenetic comparative tools for Phenome-Genome studies:
convergent amino acid substitution (CAAS) discovery, ancestral-state
scoring, relative evolutionary rate (RERconverge), directional selection
(HyPhy FADE), composite scoring, variant annotation, position-level
enrichment (POSENRICH) and functional/network enrichment.

General usage:            > nextflow run main.nf [options]
Module-specific help:     > nextflow run main.nf --help --<module_flag>

General options:
--remove_work             <true|false>  false  (the generated run script removes its work, run and
                                         Nextflow home directories after a successful run; -resume
                                         is not available afterwards)
--toy_mode, --toy_n       Smoke test on a random subset of alignments (see --help --contrast_selection)

Module flag              Description
------------------------  -------------------------------------------------
--reporting               Preliminary dataset/phenotype exploration reports
                          (independent of CAAS discovery).
--contrast_selection      Prune data and select foreground/background
                          contrasts from a continuous or discrete trait.
--ct_tool <tools>         CAAStools stages (comma-separated, e.g.
                          "discovery,resample"). resample harvests the
                          permulated labelings (the null) and can run alone;
                          discovery replays the real labeling and the null
                          over the alignments and needs those labelings.
                          See: --help --ct_tool <discovery|resample>
--ct_disambiguation       Score each CAAS by ancestral state reconstruction (ASR) and
                          label it convergent, divergent, single or no_change.
--ct_postproc             Cluster/gene-level filtering + characterization report.
--posenrich               POSENRICH: level 1 tests the share of the detections held by each
                          gene against the null cycles; the position registry, level 2 and
                          report 13 need --scoring.
--vep                     Annotate scored positions with PrimateAI-3D, COSMIC and Ensembl VEP.
                          On by default; requires --scoring (use --vep false otherwise).
--fade                    HyPhy FADE directional (top/bottom) selection scan.
--rer_tool <steps>        RERconverge relative-evolutionary-rate analysis
                          (comma-separated: build_trait,build_tree,build_matrix,continuous,binary).
--scoring                 Position- and gene-level scores from CAAS, with FADE and RER
                          evidence, calibrated against the CAAS null.
--enrichment              FCS ranked enrichment, DOMINO/STRING network analysis and
                          the comparison report. Needs --scoring; without it only the
                          RER FCS report can run.

Stages after CT can start from files of a previous run: each module help lists its
"*_from" / "*_input" options. ENRICHMENT and VEP always need --scoring in the same
invocation. See the README for all options and docs/WALKTHROUGH.md for the statistics.
'''
}

// ── CT: discovery / resample ─────────────────────────────────────────────────

def discovery_help() {
    return '''
CT Discovery — Help
=============================================
Detects Candidate Amino Acid Substitutions (CAAS) from Multiple Sequence
Alignments (MSA) against a foreground/background trait split.

The permulation core replays the real labeling (b_0) and the harvested
permulated labelings over the alignments in one pass; the discovery of the real
labeling is its b_0 slice, and the permulated ones are the null. The labelings
come from the resample stage of the same run or from --resample_from, so
discovery cannot run without them. With --caas_full_perms 0 only the real
labeling is replayed (no null).

Usage:
--alignment                <"input_dir">                null
--caas_config               <"caas_config_file">         null
--ali_format                <"fasta">                    "fasta"
--patterns                  <"1,2,3,4">                  "1,2,3"
--min_divergent_fraction    <FLOAT 0-1>                  0.5
--min_non_overlapping       <INTEGER>                    2      (minimum non-overlapping species required on FG or BG to define a CAAS)
--strict_conservation_output <true|false>   true   (with min_divergent_fraction < 1, also write the observed tables and the null of the hypotheses with no conserved pair under <outdir>/strict/)
--max_bg_gaps_fraction      <FLOAT 0-1>                  0.0
--max_fg_gaps_fraction      <FLOAT 0-1>                  0.0
--max_gaps_fraction         <FLOAT 0-1>                  0.0
--max_bg_miss_fraction      <FLOAT 0-1>                  0.0
--max_fg_miss_fraction      <FLOAT 0-1>                  0.0
--max_miss_fraction         <FLOAT 0-1>                  0.0
--miss_pair                 <true|false>                 true
--ct_core_batch_size        <INTEGER>                    20    (permulation core: genes per replay + ASR replay task; 1 = one task per gene)
--caas_export_single       <BOOL>                       false (export global caas_perms.rds for compatibility; blocks + roster are canonical)
--caas_cycle_block_size    <INTEGER>                    500   (permulation core: cycles per block; a task replays its genes against one block, the observed labeling against none; one merge task per block)
--gene_stat                max                         max            (only supported gene statistic: largest positive position score)
--publish_intermediates     <true|false>                 false

NOTE: fractions are resolved to floor(n_pairs * fraction) at runtime.
'''
}

def resample_help() {
    return '''
CT Resample — Help
=============================================
Harvests the permulated foreground/background labelings that make the CAAS null.
The trait is simulated under BM (or OU) on the species tree, its ranks are mapped
back onto the observed values, and each labeling is kept if its contrast pairs
are Dunn-independent. It reads no alignments and can run alone, for example to
inspect the null before replaying it with the discovery stage.

Usage:
--tree                      <"nwtree_file">              null
--traitname                  <"column_name">               null   (trait column of --my_traits)
--perm_strategy              <"auto|OU|BM">               "auto"  (evolutionary model of the simulation)
--perm_match_pss             <true|false>                 true   (each permulation builds its canonical pairs so that their PSS follows the PSS of the observed pairs; accepted only if every pair is within the tolerance. Not applied to count traits)
--perm_match_pss_tol          <FLOAT>                      0.25   (relative PSS tolerance of that matching)
--resample_use_n             <true|false>                 true   (use N/C sample size count traits for Jeffreys CI filtering if available)
--max_fop                    <INTEGER>                    100    (alternative hypotheses per contrast, observed and null)
--chunk_size                  <INTEGER>                    500

Permulation sizing:
--max_tries                  <INTEGER>                    1000000 (draw budget; raised 50% up to twice if the pool is short, then fails)
--caas_full_perms            <INTEGER>                    1000    (accepted permulations harvested AND replayed as the null; 0 = replay only the real labeling)

Output: directory of resample_*.tab files (one per chunk_size cycles),
permulation_manifest.tsv (tier, Dunn index and design of the canonical pairs of
each cycle) and permulation_harvest.tsv (draws tried, rejections by reason,
acceptance by PSS tolerance, pool and capacity of a sample of the draws).
'''
}

// ── Contrast selection ──────────────────────────────────────────────────────

def contrast_selection_help() {
    return '''
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
--pair_center_gate               <"none|median|mean">          "none" (candidate pairs must straddle the trait median or mean; observed selection and permulation null)
--perm_strategy                  <"auto|OU|BM">               "auto"  (evolutionary model for PSS / permulation)
--max_contrasts                   <INTEGER>                     0     (0 = dynamic discovery)
--min_contrasts                    <INTEGER>                     3     (min foreground pairs required)
--prune_data                        <true|false>                   false  (requires --reporting)
--prune_list                         <"species_list">               null
--prune_list_secondary                 <"species_list">               null
--toy_mode                           <true|false>                    false  (smoke test on a random subset of alignments; species are never subsampled)
--toy_n                               <INTEGER>                      1000   (alignments kept in --toy_mode: the first N of a shuffle seeded by --seed; an empty value or 0 falls back to 50)
'''
}

// ── CT disambiguation / post-processing ─────────────────────

def disambiguation_help() {
    return '''
CT Disambiguation — Help
=============================================
Scores each discovered CAAS with ancestral state reconstruction (ASR): how likely the
shared residue arose independently in the contrast domains. Positions are labeled
convergent, divergent, single or no_change.

Usage:
--meta_caas_from                <"meta_caas_output">          null  (a previous run's meta_caas.tsv: locates that run's outdir to reuse its permulation exports)
--ct_disambig_asr_model             <"lg">                          "lg"
--ct_disambig_asr_cache_dir          <"cache_dir">                   null
--ct_disambig_posterior_threshold      <FLOAT 0-1>                     0.1
--ct_disambig_max_tasks_per_child       <INTEGER>                       50
--asr_diagnostics                        <true|false>                    true   (render the ASR diagnostics report, 9)
'''
}

def postproc_help() {
    return '''
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
--filter_min_positions        <INTEGER>                      null (minimum candidates per interval; unset preserves legacy count gate per gene)
--filter_maxcaas                          <FLOAT 0-1>                          0.7   (filter mode)
--caas_map_dir                              <"map_dir">                          ""    (trimmer MAP tables: trains in untrimmed columns, observed and null; required by --vep)
--remove_caas_clusters                      <BOOLEAN>                            true  (discard spatial clusters)
--gene_filter_mode                          <"none|stratified">                   "stratified"
--gene_filter_strata                         <INTEGER>                            10    (gene-length strata)
--gene_filter_k_far                           <FLOAT>                              3.0   (far fence: a gene above it is removed)
--gene_filter_k_near                          <FLOAT>                              1.5   (near fence: a gene above it is removed only with a cluster train)
'''
}

// ── POSENRICH ────────────────────────────────────────────────────────────────

def posenrich_help() {
    return '''
POSENRICH — Help
=============================================
Position-level enrichment of CAAS detections on the universe U of positions (the same U as
SCORING). The universe table (universe.tsv.gz) is shared by both levels. Level 1 (Accumulation)
is an exact hypergeometric test of whether a gene holds more CAAS positions than its share of U
predicts; it runs once per direction (top/bottom/all). The position registry, its characterization, level 2 (concentration of the
signal inside genes, calibrated with the null cycles) and the overlap and proximity components read SCORING's
tables and need --scoring. Without --scoring only posenrich/universe and posenrich/accumulation are written; with
it the run adds registry, characterization, within_gene and proximity, and report 13
(13.POSENRICH_report_<trait>.html). Level 1 needs --ct_postproc or --accumulation_background_input.

Usage:
--posenrich                          <true|false>                   false
--accumulation_fdr                   <FLOAT 0-1>                    0.1   (q threshold marked in report 13 for levels 1 and 2; q itself does not depend on it)
--accumulation_caas_input            <"filtered_discovery.tsv">     null  (standalone entry: observed positives of U)
--accumulation_background_input      <"background_file">            null  (standalone entry: gene list of U)
--posenrich_background_file          <"background_file">            null  (tested positions, background.output)
--posenrich_domains, --posenrich_ucr <true|false>                   true  (blank tables are generated on the productive alignments)
--domain_variability_file, --ucr_positions_file   (optional precomputed annotation sources; each needs an explicit coordinate_space)
--pfam_cache_dir                     <"pfam_cache_dir">             ~/.cache/phylophere/pfam
--domain_ref_species                 <"species">                    ref_species_name
--fubar_sites_file                   (optional; cannot be generated in-house
                                      -- see github.com/nozerorma/ortholog_characterizator)
--posenrich_custom_file              (custom position sets; TSV: gene, zero-based CAAS position, term, description; no header)
--posenrich_collections              <"id1,id2">                    "" (level 2 collections; blank = those informative inside a protein)
--posenrich_min_genes                <int>                          2 (genes of U a term needs members in to be tested by level 2)
--posenrich_holdout                  <int>                          0 (validation: null cycles played as the observed one; 0 = off)
--posenrich_within_gene_memory       <"16 GB">                      "" (level 2 memory; blank = 16 GB)
--posenrich_radios                   <LIST>                         0,5,10 (descriptive protein-proximity radii)
--position_registry_memory, --position_module_memory   (blank sizes them from the inputs)
'''
}

// ── VEP ──────────────────────────────────────────────────────────────────────

def vep_help() {
    return '''
VEP Characterization — Help
=============================================
Annotates the positions of SCORING with PrimateAI-3D pathogenicity scores,
COSMIC somatic-mutation evidence, and/or Ensembl VEP consequence prediction
(each source skipped independently if its param is empty/unset).
Requires --scoring. Annotations are descriptive and do not change any score.

VEP reads the per-gene MAP directory from --caas_map_dir (post-processing parameter).
It has no in-house generation path: it maps each alignment codon
column to its real hg38 genomic/protein coordinate, which needs the full
alignment-to-protein pipeline, not just the alignment plus a public DB.
Generate it with: https://github.com/nozerorma/ortholog_characterizator

Usage:
--vep_primateai_db                 <"primateai_db">          ""   (PrimateAI-3D hg38 table; not bundled)
--cosmic_db                          <"cosmic_db">              ""   (COSMIC Mutant Census GRCh38 table; not bundled)
--vep_ensembl                         <BOOLEAN>                 false (runs the official Ensembl VEP CLI, independent of the two DBs above; needs --gene_ensembl_file)
--vep_mode                             <"mode">                   "auto"  (auto = database when vep_db_host is reachable, else offline cache; database; offline)
--vep_db_host                           <"host">                   "ensembldb.ensembl.org"
--vep_db_port                            <PORT>                    3306
--vep_cache_dir                        <"vep_cache_dir">          ""    (offline mode; blank auto-populates ~/.cache/phylophere/vep/<species>_<assembly> via vep_install on first use)
--vep_species                            <"species">                "homo_sapiens"
--vep_assembly                             <"assembly">               "GRCh38"
'''
}

// ── FADE ─────────────────────────────────────────────────────────────────────

def fade_help() {
    return '''
FADE (directional selection) — Help
=============================================
Runs HyPhy FADE to test directional substitution bias towards amino acids
on foreground branches in a rooted tree. Technical failures stop complete output.

Usage:
--fade_direction                         <"top"|"bottom"|"both">   "both"  (which foreground direction(s) to run)
--fade_background_scope                  <"all"|"opposite">        "all"   ("opposite" prunes the tree+alignment to foreground + opposite-extreme species only, excluding "middle" species entirely; "all" is HyPhy's own implicit background-is-everything-else behavior. With --fade_direction both, "opposite" runs each extreme group against the other.)
--fade_internal_nodes                    <"all_descendants"|"none"> "all_descendants" ("all_descendants" labels every internal branch whose descendant tips are all foreground, so all branches of a foreground clade, its stem included, are tested (HyPhy LabelTrees default); "none" tests foreground terminal branches only)
--fade_species_file                      <"path">                  null    (optional user-supplied fg/bg species file, candidate_species.tab format: species<TAB>contrast_group[<TAB>pair]; contrast_group==1 -> top, ==0 -> bottom; bypasses contrast selection for FADE when set. Requires --tree when used standalone.)
--selection_prep_batch_size          <INTEGER>                500
--fade_batch_size                      <INTEGER>                200
--fade_bf_threshold                      <FLOAT>                100  (Bayes Factor cutoff)
--fade_model                               <"LG|WAG|JTT|...">       "LG" (HyPhy built-in model; Dayhoff is unsupported)
--fade_method                                <"Variational-Bayes|Collapsed-Gibbs|Metropolis-Hastings">  "Variational-Bayes"
--fade_grid                                    <INTEGER>                20
--fade_chains, --fade_chain_length,
--fade_burn_in, --fade_samples               (MCMC-only settings; ignored by Variational-Bayes)
--fade_concentration                            <FLOAT>                  0.5
--fade_min_genes_for_heatmap                       <INTEGER>                2     (legacy compatibility option; no effect)
--fade_universe_file                              <"gene_list">            null  (else falls back to cleaned background)
--fade_tables_dir_top / --fade_tables_dir_bottom    <"dir/selection/fade/<direction>">   null  (scoring and report from the merged FADE tables of a prior run)
'''
}

// ── RERconverge ──────────────────────────────────────────────────────────────

def rer_help() {
    return '''
RERconverge (RER) — Help
=============================================
Computes Relative Evolutionary Rate — branch-length deviation correlated with
the trait — auto-routing to continuous or binary correlation.

Usage:
--rer_tool                    <"build_trait,build_tree,build_matrix,continuous,binary">   ""
--rer_trait_mode                 <"auto|continuous|binary">                          "auto"
--rer_minsp                        <INTEGER>                                            15
--winsorizeRER, --winsorizeTrait     <INTEGER>                                            3
--rer_perm_batches                    <INTEGER>                                            10   (0 = skip permulation)
--rer_perms_per_batch                   <INTEGER>                                            100
--rer_batch_seed                        <INTEGER>                                            1998 (batch i: seed + i - 1; fixed boundaries required for reproducibility)
--rer_binary_clade                        <"all|ancestral|terminal">                           "all"
--rer_min_pos                               <INTEGER>                                            2
--rer_transform                                 <"auto|ha_logit|logit|arcsin|log10|none">            "ha_logit"
--rer_universe_file                               <"gene_list">                                        null  (else falls back to cleaned background)
--rer_continuous_file                               <"trait.continuous.output">                          null  (validated tables + report from one precomputed RDS)
--rer_binary_file                                    <"trait.binary.output">                              null  (same, for binary traits; mutually exclusive with continuous file)
--rer_perms_file                                      <"perms.rds|corStat block dir">                  null  (standalone null)
--rer_gene_scores                                       <"fcs_stats.tsv">                                   null  (SCORING gene table joined to the RER gene lists)
--rer_top_n_labels                                       <INTEGER>                                            20    (genes labelled in plots)
--gene_trees, --trait_out, --trees_out, --matrix_out   (RER input and intermediate objects; see conf/rerconverge.config)
'''
}

// ── Scoring ──────────────────────────────────────────────────────────────────

def scoring_help() {
    return """
CAAS Scoring — Help
=============================================
Computes CAAS position and gene scores from CT_POSTPROC, calibrates them against the
CAAS permulation null (empirical p, Benjamini-Hochberg q) and joins FADE and RER
evidence. VEP runs after SCORING and annotates its positions.

Usage:
--scoring                   <true|false>            false
--scoring_ami                <true|false>            true    (downstream DOMINO/AMI active module analysis)
--scoring_string             <true|false>            true    (read only when --scoring_ami is unset)
--gene_stat                  <"max">                 max     (only supported gene statistic: largest positive position score)
--auto_generate_ensembl      <true|false>            false   (build --gene_ensembl_file from BioMart; read by bin/resolve_core_inputs.py)
--ensembl_dataset            <"dataset">             ""      (BioMart dataset; blank derives it from --ref_species_name)
--scoring_compare_fdr          <FLOAT 0-1>             0.15
--scoring_compare_top_n           <INTEGER>               20
--evidence_p_threshold <FLOAT 0-1> 0.05 (nominal empirical activity)
--evidence_fdr_threshold <FLOAT 0-1> 0.05 (BH highlighting)
--scoring_window_size_bp                       <INTEGER>               1000000
--caas_evidence_top_n                          <INTEGER>               0      (evidence table of the N best positions of position_scores.tsv, written after SCORING to scoring/evidence/: MRCA, tip species and residues and the score of each domain of each hypothesis; 0 = off; needs --scoring, the observed discovery.tab and ct_disambig_asr_cache_dir)
--caas_evidence_discovery                   <"discovery.tab">       ""     (discovery.tab for the evidence table when the observed tables are given and not recomputed; set by the generated reuse scripts)

Standalone mode (provide inputs directly):
--scoring_postproc_input        <"filtered_discovery.tsv">        ""
--scoring_background_input        <"cleaned_background_main.txt">   ""
--scoring_hypotheses_pairs          <"contrast_hypotheses_pairs.tsv">  ""    (domain weights; found automatically in outdir/data_exploration/2.CT/1.Traitfiles/)
--scoring_fade_summary_top          <"fade_summary_top.tsv">           ""
--scoring_fade_summary_bottom          <"fade_summary_bottom.tsv">        ""
--scoring_fade_site_top                  <"fade_site_bf_top.tsv">           ""
--scoring_fade_site_bottom                  <"fade_site_bf_bottom.tsv">        ""
--scoring_rer_input                           <"rerconverge_summary_<trait>.tsv">  ""
--scoring_vep_primateai                            <"primateai_mapped.tsv">            ""    (also the precomputed-VEP fallback of the comparison report, see enrichment.nf)
--scoring_vep_cosmic                                  <"cosmic_scores.tsv">                ""    (also the precomputed-VEP fallback of the comparison report, see enrichment.nf)
--scoring_vep_ensembl                                <"ensembl_vep_mapped.tsv">          ""    (also the precomputed-VEP fallback of the comparison report, see enrichment.nf)
--caas_null_dir                                          <"caas_permulation/">          ""    (caas_permulation/ directory of a prior run: the null of SCORING and POSENRICH)
--caas_perm_harvest_file, --caas_perm_manifest_file   <"permulation_harvest.tsv", "permulation_manifest.tsv">  ""  (null harvest and design panels of report 10)

Main outputs (scoring/): position_scores.tsv, gene_scores.tsv, fcs_stats*.tsv, gene_lists/,
position_lists/ and report 10. Gene columns: gene_caas_score (+ _top, _bottom), caas_p.gen,
rer_min_pval, rer_acceleration, fade_max_bf_<direction>.
"""
}

// ── Enrichment ───────────────────────────────────────────────────────────────

def enrichment_help() {
    return '''
Enrichment (FCS / DOMINO-STRING) — Help
=============================================
Ranked gene-set enrichment (FCS, Wilcoxon-AUC on GMT gene sets) for the CAAS and RER
rankings, DOMINO active-module detection with STRING labels, and a cross-module
Comparison report that reads the POSENRICH position characterization when it exists.
Needs --scoring; without it only the RER FCS report can run (from RER output).

Usage:
--enrichment                        <true|false>          false
--fcs_enabled                         <true|false>             true   (CAAS and RER FCS reports; does not gate AMI)
--gmt_dir                            <"gmt_dir">              <projectDir>/subworkflows/ENRICHMENT/dat
--auto_fetch_gmt                      <true|false>             false  (GUI-generated scripts only: add current GO and WikiPathways sets)
--fcs_min_genes                        <INTEGER>                10
--fcs_max_genes                        <INTEGER>                0     (0 = no cap)
--fcs_fdr                                <FLOAT 0-1>              0.15   (BH on the permulation p of each set, per ranking and database: q.perm)
--fcs_pperm_thr                            <FLOAT 0-1>              0.025  (raw permulation-p threshold of the Redlich et al. rule)
--fcs_analytic_levels                        <LIST>                   0.05,0.1,0.15  (analytic rule of Redlich et al., column redlich_level: BH of the analytic p below each level and p.perm below --fcs_pperm_thr)
--fcs_cluster_threshold                     <FLOAT 0-1>              0.2    (minimum Jaccard of the full gene sets inside a term block; complete linkage)
--fcs_top_n                                  <INTEGER>                20
--fcs_batch_size                             <INTEGER>                4 (GMT databases per task)
--fcs_block_batch_genes                      <INTEGER>                1000 (GENE_P genes per task; all cycles retained for these genes)
--caas_permulation_enrichment                    <true|false>             true

DOMINO / STRING:
--string_species                    <NCBI taxon id>          9606   (taxon of the report 12 STRING labels; the DOMINO network is human)
--string_db_dir                           <"dir">                  ""     (pre-downloaded STRING files, checked before any download)
--string_cache_dir                          <"dir">                  ""     (download cache; blank = ~/.cache/phylophere/string)
--ami_description_fdr                         <FLOAT 0-1>              0.10   (FDR of the descriptive STRING labels)
--ami_max_blocks                                <INTEGER>                8      (blocks shown in labels and summary)
--ami_block_threshold                             <FLOAT 0-1>              0.2    (minimum Jaccard to merge terms into a block)
--publish_domino_intermediates                      <true|false>             false
--domino_network_score_thr                  <INTEGER>                700
--domino_slice_thr                             <FLOAT 0-1>              0.3
--domino_module_thr                              <FLOAT 0-1>              0.05
'''
}

workflow HELP {
    // Check if --help is provided
    if (params.help) {
        // Most module toggles default to `true`/non-empty in conf/common.config and
        // conf/ct.config (the "core" CT+disambiguation+postproc+VEP chain
        // runs unless explicitly disabled), so there is no boolean signal that reliably
        // means "the user explicitly asked for this module." Instead: always show the
        // general overview first, then append detail for every module that WOULD run
        // given the current params — i.e. help for the run you are about to launch.
        // To see only the overview, explicitly disable the modules you don't need
        // (e.g. --ct_tool "" --ct_disambiguation false --ct_postproc false
        // --posenrich false --vep false --help).
        log.info general_help()

        if (params.contrast_selection)   log.info contrast_selection_help()

        // --ct_tool accepts a comma-separated list (e.g. "discovery,resample")
        if (params.ct_tool) {
            def tools = params.ct_tool.toString().split(',').collect { it.trim() }
            if ('discovery' in tools) log.info discovery_help()
            if ('resample'  in tools) log.info resample_help()
        }

        if (params.ct_disambiguation)    log.info disambiguation_help()
        if (params.ct_postproc)          log.info postproc_help()
        if (params.posenrich)            log.info posenrich_help()
        if (params.vep)                  log.info vep_help()
        if (params.fade)                 log.info fade_help()
        if (params.rer_tool)             log.info rer_help()
        if (params.scoring)              log.info scoring_help()
        if (params.enrichment)           log.info enrichment_help()

        exit 1
    }
}
