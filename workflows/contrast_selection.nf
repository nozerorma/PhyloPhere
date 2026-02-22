#!/usr/bin/env nextflow

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
# File: reporting.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CONTRAST_SELECTION Workflow: Preliminary reporting pipeline for trait analysis Rmarkdowns.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { DATASET_EXPLORATION } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_dataset_exploration"
include { DATASET_PRUNE } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_data_prune"
include { REPORTING } from "${baseDir}/workflows/reporting"
include { CI_COMPOSITION_REPORT } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ct_ci"
include { CONTRAST_ALGORITHM } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ct_independent-contrasts"
include { CHECK_MIN_CONTRASTS } from "${baseDir}/subworkflows/CT/ct_check_min_contrasts"

workflow CONTRAST_SELECTION {
    assert params.my_traits : "Contrast selection workflow requires --traitfile."
    assert params.tree : "Contrast selection workflow requires --tree."

    def trait_file = file(params.my_traits)
    def tree_file = file(params.tree)

    def dataset_out
    def reporting_out = null

    if (params.reporting && params.prune_data) {
        log.info "Reporting + pruning enabled; running a single prune/exploration pass for contrast selection."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out.pruned_results_dir)
        dataset_out = dataset_exploration_out.results_dir
    } else if (params.reporting){
        log.info "stats_df generated during reporting"
        reporting_out = REPORTING()
        dataset_out = reporting_out.dataset_out
    } else if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before contrast selection."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out.pruned_results_dir)
        dataset_out = dataset_exploration_out.results_dir
    } else {
        log.info "No stats_df provided. Rerunning dataset exploration for stats generation."
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, file('NO_FILE'))
        dataset_out = dataset_exploration_out.results_dir
    }

    // Always run composition step — the Rmd branches internally:
    // Jeffreys CI when n_trait/c_trait columns exist in the data, discrete categorization otherwise.
    log.info "Running composition analysis (CI or discrete). n_trait=${params.n_trait ?: '<none>'}, c_trait=${params.c_trait ?: '<none>'}"
    ci_composition_out = CI_COMPOSITION_REPORT(trait_file, tree_file, dataset_out)
    ci_out = ci_composition_out.results_dir

    contrast_out = CONTRAST_ALGORITHM(trait_file, tree_file, ci_out)

    // Gate: skip the entire CT pipeline if the traitfile contains fewer than
    // params.min_contrasts (default 3) foreground pairs.  When the threshold
    // is not met, CHECK_MIN_CONTRASTS writes a low_contrasts.skip sentinel to
    // outdir and emits nothing — all downstream processes are silently skipped.
    check_out = CHECK_MIN_CONTRASTS(
        contrast_out.trait_file_out,
        contrast_out.bootstrap_trait_file_out
    )

    emit:
        trait_file_out           = check_out.traitfile_out
        bootstrap_trait_file_out = check_out.boot_traitfile_out
        tree_file_out            = contrast_out.tree_file_out
        contrast_results_dir     = contrast_out.contrast_results_dir
        low_contrasts_skip       = check_out.skip_flag
}
