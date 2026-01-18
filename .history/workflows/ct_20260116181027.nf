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
#                                      
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: ct.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT Workflow: This workflow integrates the discovery and resampling modules for CAAStools.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { DISCOVERY } from "${baseDir}/subworkflows/CT/ct_discovery"
include { RESAMPLE } from "${baseDir}/subworkflows/CT/ct_resample"
include { BOOTSTRAP } from "${baseDir}/subworkflows/CT/ct_bootstrap"


// Main workflow

workflow CT {
    if (params.ct_tool) {

        def toolsToRun = params.ct_tool.split(',')

        // Initialize variables
        def discovery_out = Channel.empty()
        def resample_out = params.resample_out ?: null
        def bootstrap_out

        if (toolsToRun.contains('discovery')) {
            // Define the alignment channel align_tuple
            align_tuple = Channel
                    .fromPath("${params.alignment}/*") // Recursively search all subdirectories
                    .filter { it.isFile() } // Filter out directories
                    .map { file -> tuple(file.baseName, file) }

            discovery_out = DISCOVERY(align_tuple)
        }
        if (toolsToRun.contains('resample')) {
            // Define the tree file channel
            nw_tree = file(params.tree)
            trait_val = file(params.traitvalues)

            resample_out = RESAMPLE(nw_tree, trait_val)
        }
        if (toolsToRun.contains('bootstrap')) {
            align_tuple = Channel
                    .fromPath("${params.alignment}/*") // Recursively search all subdirectories
                    .filter { it.isFile() } // Filter out directories
                    .map { file -> tuple(file.baseName, file) }

            // Pass discovery_out if discovery was run, otherwise use empty channel
            bootstrap_out = BOOTSTRAP(align_tuple, resample_out, discovery_out)
        }
    }
}
