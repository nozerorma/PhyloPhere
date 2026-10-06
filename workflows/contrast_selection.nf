#!/usr/bin/env nextflow
// contrast_selection.nf — Select the independent foreground/background contrast pairs of a trait.
// PhyloPhere | workflows/

/*
##
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: contrast_selection.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CONTRAST_SELECTION: prepares the species tree and trait table, runs the trait
 *  reports and selects the independent contrast pairs that define the foreground
 *  and background species of the CAAS discovery.
 *
 *  Steps: tree-tip curation against the alignments (NAME_CURATION, when
 *  params.ali_sp_names or params.alignment is set); optional data pruning
 *  (params.prune_data) and exploration reports (params.reporting), or a bare
 *  dataset exploration that generates the trait statistics; the composition
 *  report (CI_COMPOSITION_REPORT) and the selection of the contrast pairs
 *  (CONTRAST_ALGORITHM); and the minimum-contrast gate (CHECK_MIN_CONTRASTS).
 *
 *  Consumes:  params.my_traits (trait table), params.tree (species tree)
 *  Produces:  traitfile, permulation traitfile and traitfile directory (gated by
 *             CHECK_MIN_CONTRASTS), tree, trait statistics, candidate species,
 *             contrast results directory, the low_contrasts skip flag and, when
 *             pruning ran, the pruned trait and tree files
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// ── Includes ─────────────────────────────────────────────────────────────────

include { DATASET_EXPLORATION } from '../subworkflows/TRAIT_ANALYSIS/ta_dataset_exploration'
include { PHENOTYPE_EXPLORATION } from '../subworkflows/TRAIT_ANALYSIS/ta_phenotype_exploration'
include { DATASET_PRUNE } from '../subworkflows/TRAIT_ANALYSIS/ta_data_prune'
include { REPORTING } from './reporting'
include { CI_COMPOSITION_REPORT } from '../subworkflows/TRAIT_ANALYSIS/ct_ci'
include { CONTRAST_ALGORITHM } from '../subworkflows/TRAIT_ANALYSIS/ct_independent-contrasts'
include { CHECK_MIN_CONTRASTS } from '../subworkflows/CT/ct_check_min_contrasts'
include { NAME_CURATION } from '../subworkflows/TRAIT_ANALYSIS/ta_name_curation'

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow CONTRAST_SELECTION {
    main:
    assert params.my_traits : "Contrast selection workflow requires --my_traits."
    assert params.tree : "Contrast selection workflow requires --tree."

    def trait_file = file(params.my_traits)
    def tree_file_ch = Channel.value(file(params.tree))

    // NAME_CURATION renames the tree tips to the alignment species names and prunes the
    // tips with no alignment counterpart. It runs here whenever alignment names are
    // available, regardless of params.reporting and params.prune_data. REPORTING() curates
    // the tree too, but only for its own reports (its emit: block does not return the
    // curated tree), so without this call the contrast pairs would be selected on the raw
    // params.tree and could include a species that has no alignment coverage.
    if (params.ali_sp_names || params.alignment) {
        def tax_id_param = params.tax_id
        def tax_id_ch = tax_id_param
            ? Channel.value(file(tax_id_param))
            : Channel.value(file('NO_FILE'))
        name_curation_out = NAME_CURATION(tree_file_ch, tax_id_ch)
        tree_file_ch = name_curation_out.curated_tree
        log.info "[CONTRAST_SELECTION] NAME_CURATION enabled — using curated tree as canonical tree."
    }

    def tree_file = tree_file_ch

    def dataset_out
    def reporting_out = null
    def contrast_stats_file
    def pruned_trait_emit = Channel.empty()
    def pruned_tree_emit = Channel.empty()

    if (params.reporting && params.prune_data) {
        log.info "Reporting + pruning enabled; running a single prune/exploration pass for contrast selection."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        
        def orig_trait_file = trait_file
        def orig_tree_file = tree_file

        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        pruned_trait_emit = prune_out.pruned_trait_file
        pruned_tree_emit = prune_out.pruned_tree_file
        dataset_exploration_out = DATASET_EXPLORATION(orig_trait_file, orig_tree_file, prune_out.pruned_results_dir)
        phenotype_out = PHENOTYPE_EXPLORATION(orig_trait_file, orig_tree_file, dataset_exploration_out.results_dir)
        dataset_out = phenotype_out.results_dir
        contrast_stats_file = dataset_exploration_out.stats_file
    } else if (params.reporting){
        log.info "stats_df generated during reporting"
        reporting_out = REPORTING()
        dataset_out = reporting_out.dataset_out
        contrast_stats_file = reporting_out.stats_file
    } else if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before contrast selection."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        
        def orig_trait_file = trait_file
        def orig_tree_file = tree_file

        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        pruned_trait_emit = prune_out.pruned_trait_file
        pruned_tree_emit = prune_out.pruned_tree_file
        dataset_exploration_out = DATASET_EXPLORATION(orig_trait_file, orig_tree_file, prune_out.pruned_results_dir)
        dataset_out = dataset_exploration_out.results_dir
        contrast_stats_file = dataset_exploration_out.stats_file
    } else {
        log.info "No stats_df provided. Rerunning dataset exploration for stats generation."
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, file('NO_FILE'))
        dataset_out = dataset_exploration_out.results_dir
        contrast_stats_file = dataset_exploration_out.stats_file
    }

    // The composition report always runs; it chooses its method from the data (see
    // CI_COMPOSITION_REPORT: Jeffreys intervals for count columns, coded levels for an
    // ordinal trait, Phylogenetic Shift Score for a continuous one).
    log.info "Running composition analysis (CI or discrete). n_trait=${params.n_trait ?: '<none>'}, c_trait=${params.c_trait ?: '<none>'}"
    ci_composition_out = CI_COMPOSITION_REPORT(trait_file, tree_file, dataset_out)
    ci_out = ci_composition_out.results_dir

    contrast_out = CONTRAST_ALGORITHM(trait_file, tree_file, ci_out)

    // Gate: with fewer than params.min_contrasts (3 when unset) foreground species in the
    // traitfile, CHECK_MIN_CONTRASTS writes a low_contrasts.skip sentinel to outdir and
    // emits no traitfile, so the processes that consume it do not run.
    check_out = CHECK_MIN_CONTRASTS(
        contrast_out.trait_file_out,
        contrast_out.permulation_trait_file_out,
        contrast_out.trait_dir_out
    )

    emit:
        trait_file_out             = check_out.traitfile_out
        permulation_trait_file_out = check_out.permulation_traitfile_out
        trait_dir_out              = check_out.trait_dir_out
        tree_file_out            = contrast_out.tree_file_out
        stats_file_out           = contrast_stats_file
        // Candidate fg/bg species pool before the Dunn gate, used by FADE (written by 3.CI-composition.Rmd).
        candidate_species_out    = ci_composition_out.candidate_species_out
        contrast_results_dir     = contrast_out.contrast_results_dir
        low_contrasts_skip       = check_out.skip_flag
        pruned_trait_file        = pruned_trait_emit
        pruned_tree_file         = pruned_tree_emit
}

