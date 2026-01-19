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
        // resample_out may be a directory (new default) or a legacy single file when running bootstrap only
        def resample_out = params.resample_out ? Channel.value(file(params.resample_out)) : null
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
            // Bootstrap must run AFTER discovery completes (if discovery runs)
            if (toolsToRun.contains('discovery')) {
                // Wait for CONCAT_DISCOVERY to complete before proceeding
                CONCAT_DISCOVERY.out.tap { discovery_done }
                
                // Define the alignment channel
                align_tuple = Channel
                        .fromPath("${params.alignment}/*")
                        .filter { it.isFile() }
                        .map { file -> tuple(file.baseName, file) }

                // Join align_tuple with discovery output from published directory
                // discovery_done ensures CONCAT_DISCOVERY has completed
                align_with_discovery = align_tuple.combine(discovery_done).map { 
                    id_file, done -> 
                    def id = id_file[0]
                    def file = id_file[1]
                    def discovery_file = new File("${params.outdir}/discovery/${id}.output")
                    tuple(id, file, discovery_file.exists() ? discovery_file : [])
                }
            } else {
                // Discovery not running; use bootstrap-only mode
                // Define the alignment channel
                align_tuple = Channel
                        .fromPath("${params.alignment}/*")
                        .filter { it.isFile() }
                        .map { file -> tuple(file.baseName, file) }

                // Add empty file as the third element (no discovery file)
                align_with_discovery = align_tuple.map { id, file -> tuple(id, file, []) }
            }
            
            bootstrap_out = BOOTSTRAP(align_with_discovery, resample_out)
            
            // Concatenate all bootstrap outputs from the published directory
            CONCAT_BOOTSTRAP(
                bootstrap_out.collect().map { "${params.outdir}/bootstrap" }
            )
        }
    }
}