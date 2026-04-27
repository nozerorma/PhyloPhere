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
include { NAME_CURATION } from "${baseDir}/subworkflows/TRAIT_ANALYSIS/ta_name_curation"

workflow REPORTING {
    assert params.my_traits : "Reporting workflow requires --traitfile."
    assert params.tree : "Reporting workflow requires --tree."

    def trait_file = file(params.my_traits)
    def tree_file_ch = Channel.value(file(params.tree))

    // NAME_CURATION: normalise tree tip labels to alignment-canonical species names.
    // Runs when ali_sp_names or alignment is provided; replaces the raw tree downstream.
    if (params.ali_sp_names || params.alignment) {
        def tax_id_ch = params.tax_id
            ? Channel.value(file(params.tax_id))
            : Channel.value(file('NO_FILE'))
        name_curation_out = NAME_CURATION(tree_file_ch, tax_id_ch)
        tree_file_ch = name_curation_out.curated_tree
        log.info "[REPORTING] NAME_CURATION enabled — using curated tree as canonical tree."
    }

    def tree_file = tree_file_ch
    def reporting_stats_file
    def pruned_trait_emit = Channel.empty()
    def pruned_tree_emit = Channel.empty()

    if (params.prune_data) {
        log.info "Pruning selected; running data pruning module before reporting."
        prune_out = DATASET_PRUNE(trait_file, tree_file)
        trait_file = prune_out.pruned_trait_file
        tree_file = prune_out.pruned_tree_file
        pruned_trait_emit = prune_out.pruned_trait_file
        pruned_tree_emit = prune_out.pruned_tree_file
        dataset_exploration_out = DATASET_EXPLORATION(trait_file, tree_file, prune_out.pruned_results_dir)
        dataset_out = dataset_exploration_out.results_dir
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
