#!/usr/bin/env nextflow
// ct.nf — Prepare the permulated labelings and the replay inputs of the CAAS permulation core.
// PhyloPhere | workflows/

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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT: resamples the phenotype labelings (RESAMPLE, params.ct_tool = resample) and
 *  prepares the inputs of the permulation core, which discovers the observed CAAS as
 *  its b_0 slice (CAAS_CORE in main.nf).
 *
 *  The trait file, permulation trait file and tree come from CONTRAST_SELECTION
 *  when it ran, and from --caas_config, --my_traits and --tree otherwise. The
 *  permuted labelings come from RESAMPLE or from --resample_from.
 *
 *  Consumes:  trait file, permulation trait file, species tree (channels or params)
 *  Produces:  trait_file and tree_file (the ones used), and, when discovery or
 *             --caas_permulation_enrichment is requested, the inputs of CAAS_CORE:
 *             caas_align_tuple (alignments to replay), caas_resample_subset (b_0 and
 *             the first N permuted labelings) and caas_blocks (the labelings in blocks of
 *             cycles, with their FOP pair weights, and b_0 apart)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// ── Includes ─────────────────────────────────────────────────────────────────

include { RESAMPLE } from '../subworkflows/CT/ct_resample'
include { listAlignmentFiles; sampleAlignmentFiles } from '../subworkflows/CT/ct_alignment_files'
include { CONCAT_RESAMPLE } from '../subworkflows/CT/ct_concat'
include { CAAS_PERMS_PREP } from '../subworkflows/CT/caas_permulation'

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow CT {
    take:
        trait_file_in
        permulation_trait_file_in
        tree_file_in
        curated_trait_in   // value channel: trait table curated by NAME_CURATION, or null
    main:
        // Output channels of the emit block, defined at workflow level.
        def trait_file_emit = Channel.empty()
        // Inputs of the permulation core: the alignments to replay, the labelings subset (b_0 and the first N
        // permuted cycles) and its cycle blocks (with the FOP pair weights). They stay empty unless discovery or
        // --caas_permulation_enrichment is requested.
        def caas_align_tuple_out = Channel.empty()
        def caas_resample_subset_out = Channel.empty()
        def caas_blocks_out = Channel.empty()
        def perm_resample_dir_out = Channel.empty()
        def tree_file_emit = Channel.empty()
        
    if (params.ct_tool) {
        if (!(params.ct_tool instanceof String)) error "ct_tool requires named stages: resample,discovery"
        def toolsToRun = params.ct_tool.split(',').collect { it.trim() }.findAll { it }.unique()
        if (!toolsToRun || toolsToRun.any { !(it in ['resample', 'discovery']) })
            error "Unknown ct_tool stage; supported stages are resample and discovery."
        def enabled = { v -> v == true || v?.toString()?.toLowerCase() in ['true','1','yes'] }
        def savedNull = params.containsKey('caas_null_dir') && params.caas_null_dir
        def observedInput = params.discovery_from ||
            (params.containsKey('disambiguation_input') && params.disambiguation_input)
        def needsCore = toolsToRun.contains('discovery') ||
            (enabled(params.ct_disambiguation) && observedInput) ||
            (enabled(params.enrichment) && enabled(params.caas_permulation_enrichment) && !savedNull)
        def align_tuple = Channel.empty()
        if (needsCore) {
            def allFiles = listAlignmentFiles(params.alignment)
            if (!allFiles) error "CT replay requires a nonempty alignment directory."
            if (allFiles.collect { it.baseName }.unique().size() != allFiles.size())
                error "CT alignment files must have unique gene IDs (filename without extension)."
            if (params.toy_mode) {
                def n = (params.toy_n != null ? params.toy_n : 50) as int
                if (n < 1) error "toy_n must be positive."
                allFiles = sampleAlignmentFiles(allFiles, n, params.seed != null ? params.seed : 1998)
            }
            align_tuple = Channel.fromList(allFiles.collect { f -> tuple(f.baseName, f) })
        }
        if (trait_file_in == null && !params.caas_config) error "CT requires --caas_config or an upstream design."
        if (tree_file_in == null && !params.tree) error "CT requires a curated tree or --tree."
        trait_file_emit = trait_file_in != null ? trait_file_in : Channel.value(file(params.caas_config, checkIfExists:true))
        tree_file_emit = tree_file_in != null ? tree_file_in : Channel.value(file(params.tree, checkIfExists:true))
        def resample_dir_out = params.resample_from
            ? Channel.value(file(params.resample_from, checkIfExists:true)) : Channel.empty()
        if (toolsToRun.contains('resample')) {
            def values = permulation_trait_file_in != null ? permulation_trait_file_in : curated_trait_in
            if (values == null && params.caas_full_perms != null && params.caas_full_perms.toString() == '0')
                values = Channel.value([])
            if (values == null) {
                if (!params.my_traits) error "CT resampling requires phenotype values."
                values = Channel.value(file(params.my_traits, checkIfExists:true))
            }
            resample_dir_out = RESAMPLE(tree_file_emit, trait_file_emit, values,
                permulation_trait_file_in != null ? 'species' : (params.sp_colname ?: 'species'))
            CONCAT_RESAMPLE(resample_dir_out)
        }
        perm_resample_dir_out = resample_dir_out
        if (needsCore) {
            if (!toolsToRun.contains('resample') && !params.resample_from)
                error "CT replay requires resample or --resample_from."
            def prep = CAAS_PERMS_PREP(trait_file_emit, resample_dir_out)
            caas_align_tuple_out = align_tuple
            caas_resample_subset_out = prep.resample_subset
            caas_blocks_out = prep.blocks
        }
    }

    emit:
        trait_file = trait_file_emit
        tree_file = tree_file_emit
        // Inputs of CAAS_CORE (main.nf): the alignments to replay, the labelings subset and its cycle blocks.
        caas_align_tuple = caas_align_tuple_out
        caas_resample_subset = caas_resample_subset_out
        caas_blocks = caas_blocks_out
        // Directory of the harvest (resample_NNN.tab and permulation_manifest.tsv); empty when no harvest ran.
        resample_dir = perm_resample_dir_out
}
