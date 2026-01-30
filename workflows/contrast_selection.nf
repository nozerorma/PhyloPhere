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

workflow CONTRAST_SELECTION {
    assert params.my_traits : "Contrast selection workflow requires --traitfile."
    assert params.tree : "Contrast selection workflow requires --tree."

    def trait_file = file(params.my_traits)
    def tree_file = file(params.tree)

    def dataset_out

    if (params.reporting){
        log.info "stats_df generated during reporting"
        def reporting_out = REPORTING()
        dataset_out = reporting_out.dataset_out
    } else if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before contrast selection."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        dataset_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out.pruned_results_dir)
    } else {
        log.info "No stats_df provided. Rerunning dataset exploration for stats generation."
        dataset_out = DATASET_EXPLORATION(trait_file, tree_file, file('NO_FILE'))
    }

    if (params.n_trait && params.c_trait) {
        log.info "Running contrast selection with population data: ${params.n_trait}, cases: ${params.c_trait}. Computing CIs."
        ci_out = CI_COMPOSITION_REPORT(trait_file, tree_file, dataset_out)
    } else {
        log.info "Running contrast selection without population data. Skipping CI computation."
        ci_out = dataset_out
    }

    contrast_out = CONTRAST_ALGORITHM(trait_file, tree_file, ci_out)
    
}
