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
include { CONCAT_DISCOVERY; CONCAT_BACKGROUND; CONCAT_RESAMPLE; CONCAT_BOOTSTRAP } from "${baseDir}/subworkflows/CT/ct_concat"

// Main workflow

workflow CT {
    take:
        trait_file_in
        bootstrap_trait_file_in
        tree_file_in
    main:
        // Output channels for emit block - must be defined at workflow level
        def discovery_concat_out = Channel.empty()
        def background_concat_out = Channel.empty()
        def background_genes_out = Channel.empty()
        def bootstrap_concat_out = Channel.empty()
        
    if (params.ct_tool) {

        def toolsToRun = params.ct_tool.split(',')

        // Define the alignment channel (used by discovery and bootstrap)
        align_tuple = Channel
                .fromPath("${params.alignment}/*") // Recursively search all subdirectories
                .filter { it.isFile() } // Filter out directories
                .map { file -> tuple(file.baseName, file) }

        // Initialize variables
        def trait_file_out
        def bootstrap_trait_file_out

        def discovery_out = Channel.empty()
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

        if (params.contrast_selection && trait_file_in && bootstrap_trait_file_in) {
            log.info "Using contrast selection output for CT analyses."
            trait_file_out = trait_file_in
            trait_val = bootstrap_trait_file_in
            tree_file_out = tree_file_in
        } else {
            log.info "No contrast selection output provided for CT analyses."
            assert params.caas_config : "CT workflow requires --caas_config."
            trait_file_out = file(params.caas_config)
            if (toolsToRun.contains('resample') || toolsToRun.contains('bootstrap')) {
                if (params.traitvalues) {
                    trait_val = file(params.traitvalues)
                    tree_file_out = file(params.tree)
                }
            }
        }

        if (toolsToRun.contains('discovery')) {
            discovery_out = DISCOVERY(align_tuple, trait_file_out)
            
            // Concatenate discovery outputs - collect actual files for staging
            discovery_out.discovery_out
                .map { id, file -> file }
                .collect()
                .ifEmpty([])
                .set { discovery_files_to_concat }
            CONCAT_DISCOVERY(discovery_files_to_concat)
            discovery_concat_out = CONCAT_DISCOVERY.out.discovery_concat
            
            // Concatenate background outputs - collect actual files for staging
            discovery_out.background_out
                .collect()
                .ifEmpty([])
                .set { background_files_to_concat }
            CONCAT_BACKGROUND(background_files_to_concat)
            background_concat_out = CONCAT_BACKGROUND.out.background_concat
            background_genes_out = CONCAT_BACKGROUND.out.background_genes
        }
        if (toolsToRun.contains('resample')) {
            // Define the tree file channel
            def resample_trigger = toolsToRun.contains('discovery')
                ? discovery_out.discovery_out.ifEmpty { Channel.value(null) }.collect()
                : Channel.value(true)

            // Handle channels differently based on whether they come from contrast_selection
            if (params.contrast_selection && trait_file_in && bootstrap_trait_file_in) {
                // tree_file_out, trait_file_out, and trait_val are already channels from CONTRAST_SELECTION
                nw_tree = tree_file_out
                caas_config = trait_file_out
                trait_values = trait_val
            } else {
                // tree_file_out, trait_file_out, and trait_val are file objects that need to be channelized
                nw_tree = resample_trigger.map { tree_file_out }
                caas_config = resample_trigger.map { trait_file_out }
                trait_values = resample_trigger.map { trait_val }
            }
            resample_out = RESAMPLE(nw_tree, caas_config, trait_values)
            
            // Concatenate all resample outputs - pass directory for staging
            CONCAT_RESAMPLE(resample_out)
        }
        if (toolsToRun.contains('bootstrap')) {
            if (toolsToRun.contains('discovery')) {
                align_with_discovery = align_tuple
                        .join(discovery_out.discovery_out)
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

            bootstrap_out = BOOTSTRAP(bootstrap_in, trait_file_out)
            
            // Concatenate all bootstrap outputs - collect actual files for staging
            bootstrap_out.bootstrap_out
                .map { id, file -> file }
                .collect()
                .ifEmpty([])
                .set { bootstrap_files }
            bootstrap_concat_out = CONCAT_BOOTSTRAP(bootstrap_files)
        }
    }
    
    emit:
        discovery_file = discovery_concat_out
        background_file = background_concat_out
        background_genes = background_genes_out
        bootstrap_file = bootstrap_concat_out
}
