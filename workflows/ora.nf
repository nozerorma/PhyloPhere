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
# File: ora.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT Workflow: This workflow integrates the discovery and resampling modules for CAAStools.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// Define the alignment channel align_tuple
align_tuple = Channel
                .fromPath("${params.alignment}/**/*")
                .filter { it.isFile() } // Filter out directories
                .map { file -> tuple(file.baseName, file) }

// Define the tree file channel
nw_tree = file(params.tree)

// Import local modules/subworkflows
include { DISCOVERY } from "${baseDir}/subworkflows/CT/ct_discovery" addParams(ALIGN_TUPLE: align_tuple)
include { RESAMPLE } from "${baseDir}/subworkflows/CT/ct_resample" addParams(NW_TREE: nw_tree)
include { BOOTSTRAP } from "${baseDir}/subworkflows/CT/ct_bootstrap" addParams(ALIGN_TUPLE: align_tuple)

// Main workflow
def toolsToRun = params.ct_tool.split(',')

workflow CT {
    if (toolsToRun.contains('discovery')) {
        discovery_out = DISCOVERY(align_tuple)
    }
    if (toolsToRun.contains('resample')) {
        resample_out = RESAMPLE(nw_tree)
    }
    if (toolsToRun.contains('bootstrap')) {
        boostrap_out = BOOTSTRAP(align_tuple, resample_out)
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nYay!\n" : "Oops .. something went wrong" )
}
