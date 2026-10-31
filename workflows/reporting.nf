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
 *  REPORTING Workflow: Preliminary reporting pipeline for trait analysis Rmarkdowns, run on
 *  the trait table and species tree curated by NAME_CURATION (main.nf).
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { DATASET_EXPLORATION } from '../subworkflows/TRAIT_ANALYSIS/ta_dataset_exploration'
include { PHENOTYPE_EXPLORATION } from '../subworkflows/TRAIT_ANALYSIS/ta_phenotype_exploration'
include { PRUNE_DATASET as DATASET_PRUNE } from '../subworkflows/TRAIT_ANALYSIS/ta_data_prune'

workflow REPORTING {
    take:
        curated_trait_ch   // value channel: trait table curated by NAME_CURATION (main.nf)
        curated_tree_ch    // value channel: tree curated by NAME_CURATION (main.nf)

    main:
    def pruned_trait_emit = Channel.empty()
    def pruned_tree_emit = Channel.empty()
    def prune_results = file('NO_FILE')
    if (params.prune_data) {
        def prune_out = DATASET_PRUNE(curated_trait_ch, curated_tree_ch)
        prune_results = prune_out.pruned_results_dir
        pruned_trait_emit = prune_out.pruned_trait_file
        pruned_tree_emit = prune_out.pruned_tree_file
    }
    // Reports retain the original inputs for before/after comparison; shared
    // setup loads the explicitly staged pruned files when pruning ran.
    def dataset_exploration_out = DATASET_EXPLORATION(curated_trait_ch, curated_tree_ch, prune_results)
    def phenotype_out = PHENOTYPE_EXPLORATION(curated_trait_ch, curated_tree_ch, dataset_exploration_out.results_dir)
    def reporting_stats_file = dataset_exploration_out.stats_file

    emit:
        dataset_out = phenotype_out.results_dir
        stats_file = reporting_stats_file
        pruned_trait_file = pruned_trait_emit
        pruned_tree_file = pruned_tree_emit
}
