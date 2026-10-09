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
include { DATASET_PRUNE } from '../subworkflows/TRAIT_ANALYSIS/ta_data_prune'

workflow REPORTING {
    take:
        curated_trait_ch   // value channel: trait table curated by NAME_CURATION (main.nf)
        curated_tree_ch    // value channel: tree curated by NAME_CURATION (main.nf)

    main:
    // The names are settled by NAME_CURATION, run once by main.nf: the curated tree and trait
    // table replace --tree and --my_traits here.
    def trait_file = curated_trait_ch
    def tree_file_ch = curated_tree_ch

    def tree_file = tree_file_ch
    def reporting_stats_file
    def pruned_trait_emit = Channel.empty()
    def pruned_tree_emit = Channel.empty()

    if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before reporting."
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
        reporting_stats_file = dataset_exploration_out.stats_file
    } else {
        log.info "No data pruning selected; skipping data pruning module."
        prune_out = file('NO_FILE')
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out)
        phenotype_out = PHENOTYPE_EXPLORATION(trait_file, tree_file, dataset_exploration_out.results_dir)
        dataset_out = phenotype_out.results_dir
        reporting_stats_file = dataset_exploration_out.stats_file
    }

    emit:
        dataset_out
        stats_file = reporting_stats_file
        pruned_trait_file = pruned_trait_emit
        pruned_tree_file = pruned_tree_emit
}
