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
 *  REPORTING Workflow: Preliminary reporting pipeline for trait analysis Rmarkdowns.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { DATASET_EXPLORATION } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_dataset_exploration"
include { PHENOTYPE_EXPLORATION } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_phenotype_exploration"
include { DATASET_PRUNE } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_data_prune"

workflow REPORTING {
    assert params.my_traits : "Reporting workflow requires --traitfile."
    assert params.tree : "Reporting workflow requires --tree."

    def trait_file = file(params.my_traits)
    def tree_file = file(params.tree)

    if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before reporting."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        dataset_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out.pruned_results_dir)
    } else {
        log.info "No data pruning selected; skipping data pruning module."
        prune_out = file('NO_FILE')
        dataset_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out)
        phenotype_out = PHENOTYPE_EXPLORATION(trait_file, tree_file, dataset_out)
    }

    emit:
        dataset_out
}
