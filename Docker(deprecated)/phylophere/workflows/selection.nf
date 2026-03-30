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

        // Define the alignment channel (used by discovery and bootstrap)
        align_tuple = Channel
                .fromPath("${params.alignment}/*") // Recursively search all subdirectories
                .filter { it.isFile() } // Filter out directories
                .map { file -> tuple(file.baseName, file) }

        // Initialize variables
        def discovery_out = Channel.empty()
        def discovery_done = Channel.value(true)
        // resample_out may be a directory (new default) or a legacy single file when running bootstrap only
        def resample_out = Channel.empty()
        if (params.resample_out) {
            def resample_path = file(params.resample_out)
            if (resample_path.isDirectory()) {
                resample_out = Channel.value(file(params.resample_out, type: 'dir'))
            } else {
                resample_out = Channel.value(resample_path)
            }
        }
        def bootstrap_out

        if (toolsToRun.contains('discovery')) {
            discovery_out = DISCOVERY(align_tuple)
            discovery_done = discovery_out.collect()
            
            // Concatenate all discovery outputs from the published directory
            // Wait for all DISCOVERY processes to complete, then collect from publishDir
            CONCAT_DISCOVERY(
                discovery_done.map { "${params.outdir}/discovery" }
            )
        }
        if (toolsToRun.contains('resample')) {
            // Define the tree file channel
            nw_tree = discovery_done.map { file(params.tree) }
            trait_val = discovery_done.map { file(params.traitvalues) }

            resample_out = RESAMPLE(nw_tree, trait_val)
            
            // Concatenate all resample outputs from the published directory
            // The resample output is a directory, so we need to look inside it
            CONCAT_RESAMPLE(
                resample_out.map { "${params.outdir}/resample/${it.name}" }
            )
        }
        if (toolsToRun.contains('bootstrap')) {
            if (toolsToRun.contains('discovery')) {
                align_with_discovery = align_tuple
                        .join(discovery_out)
                        .map { row -> tuple(row[0], row[1], row[2]) }
                        .combine(discovery_done)
                        .map { row -> tuple(row[0], row[1], row[2]) }
            } else if (params.discovery_out && params.discovery_out != "none") {
                def discovery_path = file(params.discovery_out)
                if (discovery_path.isDirectory()) {
                    def discovery_files = Channel
                            .fromPath("${params.discovery_out}/*")
                            .filter { it.isFile() }
                            .map { file -> tuple(file.baseName, file) }

                    align_with_discovery = align_tuple
                            .join(discovery_files)
                            .map { row -> tuple(row[0], row[1], row[2]) }
                } else {
                    def discovery_file_ch = Channel.value(discovery_path)
                    align_with_discovery = align_tuple
                            .combine(discovery_file_ch)
                            .map { row -> tuple(row[0], row[1], row[2]) }
                }
            } else {
                align_with_discovery = align_tuple.map { id, alignmentFile -> tuple(id, alignmentFile, []) }
            }

            bootstrap_in = align_with_discovery
                    .combine(resample_out)
                    .map { row -> tuple(row[0], row[1], row[2] ?: [], row[3]) }

            bootstrap_out = BOOTSTRAP(bootstrap_in)
            
            // Concatenate all bootstrap outputs from the published directory
            def bootstrap_done = bootstrap_out.bootstrap_out.ifEmpty { Channel.value(null) }.collect()
            // If discovery was run, ensure CONCAT_BOOTSTRAP waits for discovery completion too
            def concat_trigger = toolsToRun.contains('discovery') ? 
                bootstrap_done.combine(discovery_done).map { "${params.outdir}/bootstrap" } :
                bootstrap_done.map { "${params.outdir}/bootstrap" }
            CONCAT_BOOTSTRAP(concat_trigger)
        }
    }
}
