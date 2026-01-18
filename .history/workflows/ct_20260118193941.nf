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
include { CONCAT_DISCOVERY; CONCAT_RESAMPLE; CONCAT_BOOTSTRAP } from "${baseDir}/subworkflows/CT/ct_concat"


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
            
            // Concatenate all discovery outputs from the published directory
            // Wait for all DISCOVERY processes to complete, then collect from publishDir
            CONCAT_DISCOVERY(
                discovery_out.collect().map { "${params.outdir}/discovery" }
            )
        }
        if (toolsToRun.contains('resample')) {
            // Define the tree file channel
            nw_tree = file(params.tree)
            trait_val = file(params.traitvalues)

            resample_out = RESAMPLE(nw_tree, trait_val)
            
            // Concatenate all resample outputs from the published directory
            // The resample output is a directory, so we need to look inside it
            CONCAT_RESAMPLE(
                resample_out.map { "${params.outdir}/resample/${it.name}" }
            )
        }
        if (toolsToRun.contains('bootstrap')) {
            // Define the alignment channel
            align_tuple = Channel
                    .fromPath("${params.alignment}/*")
                    .filter { it.isFile() }
                    .map { file -> tuple(file.baseName, file) }

            // If discovery was run, join with discovery output; otherwise add null
            if (toolsToRun.contains('discovery')) {
                // Join align_tuple with discovery_out by alignmentID
                align_with_discovery = align_tuple.join(discovery_out)
            } else {
                // Add null as the third element (no discovery file)
                align_with_discovery = align_tuple.map { id, file -> tuple(id, file, null) }
            }
            
            // resample_out can be a string (from params.resample_out parameter) or a channel
            // Convert channel to string path for consistency
            def resample_path = resample_out instanceof String ? resample_out : resample_out?.map { it.toString() }
            
            bootstrap_out = BOOTSTRAP(align_with_discovery, resample_path)
            
            // Concatenate all bootstrap outputs from the published directory
            CONCAT_BOOTSTRAP(
                bootstrap_out.collect().map { "${params.outdir}/bootstrap" }
            )
        }
    }
}
