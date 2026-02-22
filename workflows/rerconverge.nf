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
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: rerconverge.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Unlock the secrets of evolutionary relationships with Phylly! 🌳🔍 This Nextflow pipeline
 * packs a powerful punch, offering a comprehensive suite of phylogenetic comparative tools
 * and analyses. Dive into the world of evolutionary biology like never before and elevate
 * your research to new heights! 🚀🧬 #Phylly #EvolutionaryInsights #NextflowPipeline
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { RER_TRAIT } from "${baseDir}/subworkflows/RERCONVERGE/rer_trait"
include { RER_TREES } from "${baseDir}/subworkflows/RERCONVERGE/rer_trees"
include { RER_MATRIX } from "${baseDir}/subworkflows/RERCONVERGE/rer_matrix"
include { RER_CONT } from "${baseDir}/subworkflows/RERCONVERGE/rer_cont"
//include { RER_BIN } from "${baseDir}/subworkflows/RERCONVERGE/rer_bin" addParams(TREE_TUPLE: tree_tuple)
//include { RER_ENRICH } from "${baseDir}/subworkflows/RERCONVERGE/rer_enrich" addParams(TREE_TUPLE: tree_tuple)

// Define the gene tree channel tree_tuple
gene_trees_file = file( params.gene_trees )
// Define the trait file channel
my_traitfile = file( params.my_traits )

// Main workflow
workflow RER_MAIN {

    // Define trait_out, trees_out, and matrix_out
    trait_out = params.trait_out
    masterTrees_out = params.trees_out
    matrix_out = params.matrix_out

    if (params.rer_tool) {
        def toolsToRun = params.rer_tool.split(',')

        // Conditionally run the 'build_trait' tool
        if (toolsToRun.contains('build_trait')) {
            trait_out = RER_TRAIT(my_traitfile)
        }

        // Conditionally run the 'build_tree' tool
        if (toolsToRun.contains('build_tree')) {
            trees_out = RER_TREES(my_traitfile, gene_trees_file)
            trees_out.view { file ->
                println "Output tree file: ${file}"
            }

            prunedTrees_out = trees_out[1]  // Capturing pruned trees file
            masterTrees_out = trees_out[0]   // Capturing master tree output
        }


        // Conditionally run the 'build_matrix' tool
        if (toolsToRun.contains('build_matrix')) {
            matrix_out = RER_MATRIX(trait_out, masterTrees_out)
        }

        // Conditionally run the 'continuous' tool
        if (toolsToRun.contains('continuous')) {
            // Use outputs from the 'build_trait' tool if available, otherwise use defaults from nextflow.config
            continuous_out = RER_CONT(trait_out, masterTrees_out, matrix_out)
        }
    }
}
