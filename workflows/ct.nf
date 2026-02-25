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
        def background_raw_out = Channel.empty()
        def background_genes_out = Channel.empty()
        def bootstrap_concat_out = Channel.empty()
        def trait_file_emit = Channel.empty()
        def tree_file_emit = Channel.empty()
        
    if (params.ct_tool) {
        // Guard: params.ct_tool may be a Boolean (true) when --ct_tool is
        // passed without a value by some shells/Nextflow CLI versions.
        // Coerce to String first so .split() does not trigger a DSL2 error.
        def toolsToRun = params.ct_tool instanceof String
            ? params.ct_tool.split(',').collect { it.trim() }.findAll { it }
            : []

        // Define the alignment channel (used by discovery and bootstrap)
        align_tuple = Channel
                .fromPath("${params.alignment}/*") // Recursively search all subdirectories
                .filter { it.isFile() } // Filter out directories
                .map { file -> tuple(file.baseName, file) }

        // Initialize variables
        def trait_file_out
        def bootstrap_trait_file_out

        def discovery_out = Channel.empty()
        def discovery_done = Channel.value(true)
        // resample_dir_out  → partitioned directory passed to ct bootstrap -s
        // resample_out      → concatenated resample.tab used for reporting / emit
        def resample_out = Channel.empty()
        def resample_dir_out = Channel.empty()   // directory channel for BOOTSTRAP
        if (params.resample_out) {
            def resample_path = file(params.resample_out)
            if (resample_path.isDirectory()) {
                resample_dir_out = Channel.value(file(params.resample_out, type: 'dir'))
                resample_out     = resample_dir_out
            } else {
                // legacy single-file fallback (pre-partitioned runs)
                resample_out     = Channel.value(resample_path)
                resample_dir_out = resample_out
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
            tree_file_out = file(params.tree)
            if (toolsToRun.contains('resample') || toolsToRun.contains('bootstrap')) {
                if (params.traitvalues) {
                    trait_val = file(params.traitvalues)
                }
            }
        }

        // Normalize trait/tree output channels for downstream modules
        trait_file_emit = (params.contrast_selection && trait_file_in && bootstrap_trait_file_in) ? trait_file_out : Channel.value(trait_file_out)
        tree_file_emit  = (params.contrast_selection && trait_file_in && bootstrap_trait_file_in) ? tree_file_out  : Channel.value(tree_file_out)

        if (toolsToRun.contains('discovery')) {
            discovery_out = DISCOVERY(align_tuple, trait_file_out)

            // Hard barrier: downstream steps should start only after discovery is fully complete
            discovery_done = discovery_out.discovery_out
                .collect()
                .ifEmpty([])
            
            // Concatenate discovery outputs - collect actual files for staging
            discovery_out.discovery_out
                .map { id, file -> file }
                .collect()
                .ifEmpty([])
                .set { discovery_files_to_concat }
            CONCAT_DISCOVERY(discovery_files_to_concat)
            discovery_concat_out = CONCAT_DISCOVERY.out.discovery_concat
            
            // Concatenate background outputs - collect actual files for staging
            background_raw_out = discovery_out.background_out
            discovery_out.background_out
                .collect()
                .ifEmpty([])
                .set { background_files_to_concat }
            CONCAT_BACKGROUND(background_files_to_concat)
            background_concat_out = CONCAT_BACKGROUND.out.background_concat
            background_genes_out = CONCAT_BACKGROUND.out.background_genes
        }
        if (toolsToRun.contains('resample')) {
            // Discovery barrier trigger (if discovery was requested, wait until it fully completes)
            def resample_trigger = discovery_done

            // Handle channels differently based on whether they come from contrast_selection
            if (params.contrast_selection && trait_file_in && bootstrap_trait_file_in) {
                // tree_file_out, trait_file_out, and trait_val are already channels from CONTRAST_SELECTION.
                // Combine with resample_trigger to enforce discovery → resample ordering.
                // File-staging collisions are prevented by stageAs aliases in the RESAMPLE process.
                nw_tree = tree_file_out
                    .combine(resample_trigger)
                    .map { row ->
                        (row instanceof List || row instanceof Object[]) ? row[0] : row
                    }
                caas_config = trait_file_out
                    .combine(resample_trigger)
                    .map { row ->
                        (row instanceof List || row instanceof Object[]) ? row[0] : row
                    }
                trait_values = trait_val
                    .combine(resample_trigger)
                    .map { row ->
                        (row instanceof List || row instanceof Object[]) ? row[0] : row
                    }
            } else {
                // tree_file_out, trait_file_out, and trait_val are file objects that need to be channelized
                nw_tree = resample_trigger.map { file(tree_file_out) }
                caas_config = resample_trigger.map { file(trait_file_out) }
                trait_values = resample_trigger.map { file(trait_val) }
            }
            resample_dir_out = RESAMPLE(nw_tree, caas_config, trait_values)

            // Concatenate the partitioned directory into a single resample.tab for reporting
            CONCAT_RESAMPLE(resample_dir_out)
            resample_out = CONCAT_RESAMPLE.out.resample_concat
            // NOTE: resample_dir_out retains the raw directory so bootstrap receives
            // the partitioned resample_NNN.tab files, not the merged flat file.
        }
        if (toolsToRun.contains('bootstrap')) {
            if (toolsToRun.contains('discovery')) {
                // Use discovery results from the pipeline
                align_with_discovery = align_tuple
                        .join(discovery_out.discovery_out)
                        .map { row -> tuple(row[0], row[1], row[2]) }

                // Hard barrier: bootstrap starts only after full discovery completion
                align_with_discovery = align_with_discovery
                        .combine(discovery_done)
                        .map { row -> tuple(row[0], row[1], row[2]) }
            } else if (params.discovery_out && params.discovery_out != "none") {
                // Use external discovery file(s)
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
                    // Single concatenated discovery.tab:
                    // Read gene names at plan-time, filter align_tuple to only those
                    // genes, then pass the full discovery file to bootstrap --discovery.
                    // Alignment IDs have the form GENE.Homo_sapiens.filter2, so we
                    // extract the gene name as the token before the first dot.
                    def caas_genes = discovery_path
                        .readLines()
                        .drop(1)                          // skip header
                        .collect { it.split('\t')[0] }   // Gene column
                        .toSet()

                    align_with_discovery = align_tuple
                        .filter { id, f -> caas_genes.contains(id.split('\\.')[0]) }
                        .map    { id, f -> tuple(id, f, discovery_path) }
                }
            } else {
                // No discovery file - create placeholder for bootstrap process
                // Use a marker file to indicate no discovery optimization
                align_with_discovery = align_tuple.map { id, alignmentFile -> 
                    tuple(id, alignmentFile, file('NO_FILE')) 
                }
            }

            // Pass the partitioned resample directory (resample_dir_out), NOT the
            // concatenated resample.tab, so ct bootstrap -s receives individual
            // resample_NNN.tab files as expected.
            bootstrap_in = align_with_discovery
                    .combine(resample_dir_out)
                    .map { row -> tuple(row[0], row[1], row[2], row[3]) }

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
        background_file_raw = background_raw_out
        background_file = background_concat_out
        background_genes = background_genes_out
        bootstrap_file = bootstrap_concat_out
        trait_file = trait_file_emit
        tree_file = tree_file_emit
}
