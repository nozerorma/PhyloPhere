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
include { DISCOVERY; DISCOVERY_BATCHED } from "${baseDir}/subworkflows/CT/ct_discovery"
include { RESAMPLE } from "${baseDir}/subworkflows/CT/ct_resample"
include { CONCAT_DISCOVERY; CONCAT_BACKGROUND; CONCAT_RESAMPLE } from "${baseDir}/subworkflows/CT/ct_concat"
include { CAAS_PERMS_PREP } from "${baseDir}/subworkflows/CT/caas_permulation"

// Main workflow

workflow CT {
    take:
        trait_file_in
        permulation_trait_file_in
        tree_file_in
    main:
        def createBatchManifestText = { List<String> rows ->
            rows
                .collect { row -> row.replaceFirst(/^\s+/, '') }
                .join(System.lineSeparator()) + System.lineSeparator()
        }

        // Output channels for emit block - must be defined at workflow level
        def discovery_concat_out = Channel.empty()
        def background_concat_out = Channel.empty()
        def background_raw_out = Channel.empty()
        def background_genes_out = Channel.empty()
        def trait_file_emit = Channel.empty()
        // CAAS permulation-excess null (full-pool perm-discovery export). Populated
        // only when --caas_permulation_enrichment is set and perm-replay+resample run.
        def caas_perm_discovery_out = Channel.empty()
        def caas_resample_subset_out = Channel.empty()
        def caas_fop_pairs_out = Channel.value(file('NO_FOP_PAIRS'))
        def tree_file_emit = Channel.empty()
        
    if (params.ct_tool) {
        // Guard: params.ct_tool may be a Boolean (true) when --ct_tool is
        // passed without a value by some shells/Nextflow CLI versions.
        // Coerce to String first so .split() does not trigger a DSL2 error.
        def toolsToRun = params.ct_tool instanceof String
            ? params.ct_tool.split(',').collect { it.trim() }.findAll { it }
            : []

        // Define the alignment channel (used by discovery and the permulation-excess null).
        // Accepts either a plain directory or a .tar.gz archive.
        // For archives, members are listed at startup and extracted on-demand inside each process.
        // When toy_mode=true, a random subset of toy_n alignments is used.
        def alignParam = params.alignment as String
        def allFiles = file(alignParam).listFiles()?.findAll { it.isFile() && !it.name.matches('.*\\.txt$|.*\\.tsv$|.*\\.csv$|.*\\.log$|.*\\.map$') } ?: []
        if (params.toy_mode) {
            def n = (params.toy_n ?: 50) as int
            // Seeded (not a bare `Collections.shuffle(allFiles)`) so the same gene subset
            // is picked every run with the same --seed. Without this, -resume is close to
            // useless for toy_mode runs: DISCOVERY_BATCHED cache hits
            // depend on which specific alignment files a batch contains, and an unseeded
            // shuffle picks a genuinely different random subset on every invocation, so
            // essentially nothing from a prior run's cache ever matches a resume attempt.
            Collections.shuffle(allFiles, new Random(params.seed as long))
            allFiles = allFiles.take(n)
            log.info "[toy_mode] CT: using ${allFiles.size()} randomly sampled alignments from directory (seed=${params.seed})"
        }
        align_tuple = Channel
            .fromList(allFiles.collect { f -> tuple(f.baseName, f) })

        // Initialize variables
        def trait_file_out
        def permulation_trait_file_out
        def discovery_results = Channel.empty()
        def background_results = Channel.empty()

        def discovery_out = Channel.empty()
        def discovery_done = Channel.value(true)
        // resample_dir_out  → partitioned directory passed to the permulation-excess null
        // resample_out      → concatenated resample.tab used for reporting / emit
        def resample_out = Channel.empty()
        def resample_dir_out = Channel.empty()   // directory channel for CAAS_PERMS_PREP
        if (params.resample_from) {
            def resample_path = file(params.resample_from)
            if (resample_path.isDirectory()) {
                resample_dir_out = Channel.value(file(params.resample_from, type: 'dir'))
                resample_out     = resample_dir_out
            } else {
                // legacy single-file fallback (pre-partitioned runs)
                resample_out     = Channel.value(resample_path)
                resample_dir_out = resample_out
            }
        }
        if (params.contrast_selection && trait_file_in && permulation_trait_file_in) {
            log.info "Using contrast selection output for CT analyses."
            trait_file_out = trait_file_in
            trait_val = permulation_trait_file_in
            tree_file_out = tree_file_in
        } else {
            log.info "No contrast selection output provided for CT analyses."
            assert params.caas_config : "CT workflow requires --caas_config."
            trait_file_out = file(params.caas_config)
            tree_file_out = file(params.tree)
            if (toolsToRun.contains('resample')) {
                if (params.my_traits) {
                    trait_val = file(params.my_traits)
                }
            }
        }

        // Normalize trait/tree output channels for downstream modules
        trait_file_emit = (params.contrast_selection && trait_file_in && permulation_trait_file_in) ? trait_file_out : Channel.value(trait_file_out)
        tree_file_emit  = (params.contrast_selection && trait_file_in && permulation_trait_file_in) ? tree_file_out  : Channel.value(tree_file_out)

        if (toolsToRun.contains('discovery')) {
            def discoveryBatchSize = (params.ct_discovery_batch_size ?: 1) as int
            if (discoveryBatchSize > 1) {
                def discoveryBatchCounter = 0
                def discovery_batches = align_tuple
                    .collate(discoveryBatchSize)
                    .map { batch ->
                        def batchID = sprintf('discovery_batch_%05d', ++discoveryBatchCounter)
                        def manifestText = createBatchManifestText(
                            batch.collect { row -> "${row[0]}\t${row[1].name}" }
                        )
                        tuple(batchID, batch.size(), manifestText, batch.collect { row -> row[1] }.unique())
                    }

                discovery_out = DISCOVERY_BATCHED(discovery_batches, trait_file_out)
                discovery_results = discovery_out.discovery_out
                    .flatten()
                    .map { file -> tuple(file.baseName, file) }
                background_results = discovery_out.background_out
                    .flatten()
            } else {
                discovery_out = DISCOVERY(align_tuple, trait_file_out)
                discovery_results = discovery_out.discovery_out
                background_results = discovery_out.background_out
            }

            // Hard barrier: downstream steps should start only after discovery is fully complete
            discovery_done = discovery_results
                .collect()
                .ifEmpty([])
            
            // Concatenate discovery outputs - collect actual files for staging
            discovery_results
                .map { id, file -> file }
                .collect()
                .ifEmpty([])
                .set { discovery_files_to_concat }
            CONCAT_DISCOVERY(discovery_files_to_concat)
            discovery_concat_out = CONCAT_DISCOVERY.out.discovery_concat
            
            // Concatenate background outputs - collect actual files for staging
            background_raw_out = background_results
            background_results
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
            if (params.contrast_selection && trait_file_in && permulation_trait_file_in) {
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
            // NOTE: resample_dir_out retains the raw directory so perm-replay receives
            // the partitioned resample_NNN.tab files, not the merged flat file.
        }
        // CAAS permulation-excess null: a full-pool pass over N permuted labelings,
        // replayed through analyze_gene_disambiguation downstream. Only needs
        // align_tuple/trait_file_out/resample_dir_out, all already in scope from
        // discovery/resample above — no longer gated on the (now-removed) resample CLI tool.
        if (params.caas_permulation_enrichment) {
            def perms_prep = CAAS_PERMS_PREP(align_tuple, trait_file_out, resample_dir_out)
            caas_perm_discovery_out  = perms_prep.perm_discovery
            caas_resample_subset_out = perms_prep.resample_subset
            caas_fop_pairs_out       = perms_prep.fop_pairs
        }
    }
    
    emit:
        discovery_file = discovery_concat_out
        background_file_raw = background_raw_out
        background_file = background_concat_out
        background_genes = background_genes_out
        trait_file = trait_file_emit
        tree_file = tree_file_emit
        // CAAS permulation-excess: full-pool perm-discovery + the N-cycle resample
        // subset, consumed downstream by CAAS_PERMULATION (main.nf) → caas_perms.rds.
        caas_perm_discovery = caas_perm_discovery_out
        caas_resample_subset = caas_resample_subset_out
        caas_fop_pairs = caas_fop_pairs_out
}
