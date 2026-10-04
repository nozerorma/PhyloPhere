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
 *  CT Workflow: resamples the phenotype labelings and prepares the inputs of the permulation core,
 *  which discovers the observed CAAS as its b_0 slice (CAAS_CORE, main.nf).
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { RESAMPLE } from '../subworkflows/CT/ct_resample'
include { listAlignmentFiles; sampleAlignmentFiles } from '../subworkflows/CT/ct_alignment_files'
include { CONCAT_RESAMPLE } from '../subworkflows/CT/ct_concat'
include { CAAS_PERMS_PREP } from '../subworkflows/CT/caas_permulation'

// Main workflow

workflow CT {
    take:
        trait_file_in
        permulation_trait_file_in
        tree_file_in
    main:
        // Output channels for emit block - must be defined at workflow level
        def trait_file_emit = Channel.empty()
        // Inputs of the permulation core: the alignments to replay, the labelings subset (b_0 and the first N
        // permuted cycles) and the FOP pair weights. Populated when discovery is requested or
        // --caas_permulation_enrichment is set.
        def caas_align_tuple_out = Channel.empty()
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

        // Define the alignment channel (replayed by the permulation core).
        // params.alignment is a directory of per-gene alignment files.
        // When toy_mode=true, a seeded random subset of toy_n alignments is used (the same one for the same
        // --seed), which also keeps the batches of the core, and so -resume cache hits, stable.
        def allFiles = listAlignmentFiles(params.alignment)
        if (params.toy_mode) {
            def n = (params.toy_n ?: 50) as int
            allFiles = sampleAlignmentFiles(allFiles, n, params.seed ?: 1998)
            log.info "[toy_mode] CT: using ${allFiles.size()} randomly sampled alignments from directory (seed=${params.seed})"
        }
        align_tuple = Channel
            .fromList(allFiles.collect { f -> tuple(f.baseName, f) })

        // Initialize variables
        def trait_file_out
        def permulation_trait_file_out
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

        if (toolsToRun.contains('resample')) {
            // Handle channels differently based on whether they come from contrast_selection
            if (params.contrast_selection && trait_file_in && permulation_trait_file_in) {
                // tree_file_out, trait_file_out, and trait_val are already channels from CONTRAST_SELECTION.
                // File-staging collisions are prevented by stageAs aliases in the RESAMPLE process.
                nw_tree = tree_file_out
                caas_config = trait_file_out
                trait_values = trait_val
            } else {
                // tree_file_out, trait_file_out, and trait_val are file objects that need to be channelized
                nw_tree = Channel.value(file(tree_file_out))
                caas_config = Channel.value(file(trait_file_out))
                trait_values = Channel.value(file(trait_val))
            }
            resample_dir_out = RESAMPLE(nw_tree, caas_config, trait_values)

            // Concatenate the partitioned directory into a single resample.tab for reporting
            CONCAT_RESAMPLE(resample_dir_out)
            resample_out = CONCAT_RESAMPLE.out.resample_concat
            // NOTE: resample_dir_out retains the raw directory so perm-replay receives
            // the partitioned resample_NNN.tab files, not the merged flat file.
        }
        // Permulation core: b_0 and the first N permuted labelings are subset here; CAAS_CORE (main.nf)
        // replays them over the alignments in align_tuple. Needs trait_file_out and resample_dir_out, in scope from
        // the resample step above. Discovery is the b_0 slice of that replay, so it needs the labelings too.
        if (toolsToRun.contains('discovery') || params.caas_permulation_enrichment) {
            if (toolsToRun.contains('discovery') && !toolsToRun.contains('resample') && !params.resample_from) {
                error "ct_tool 'discovery' needs the permuted labelings of the core: add 'resample' to --ct_tool or give --resample_from"
            }
            def perms_prep = CAAS_PERMS_PREP(trait_file_out, resample_dir_out)
            caas_align_tuple_out     = align_tuple
            caas_resample_subset_out = perms_prep.resample_subset
            caas_fop_pairs_out       = perms_prep.fop_pairs
        }
    }
    
    emit:
        trait_file = trait_file_emit
        tree_file = tree_file_emit
        // CAAS permulation-excess: the alignments to replay + the N-cycle resample
        // subset, consumed downstream by CAAS_CORE (main.nf) → caas_perms.rds.
        caas_align_tuple = caas_align_tuple_out
        caas_resample_subset = caas_resample_subset_out
        caas_fop_pairs = caas_fop_pairs_out
}
