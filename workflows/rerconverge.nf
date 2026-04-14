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
# File: rerconverge.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Unlock the secrets of evolutionary relationships with Phylly! 🌳🔍 This Nextflow pipeline
 * packs a powerful punch, offering a comprehensive suite of phylogenetic comparative tools
 * and analyses. Dive into the world of evolutionary biology like never before and elevate
 * your research to new heights! 🚀🧬 #Phylly #EvolutionaryInsights #NextflowPipeline
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local modules/subworkflows
include { RER_TRAIT }          from "${baseDir}/subworkflows/RERCONVERGE/rer_trait"
include { RER_TREES }          from "${baseDir}/subworkflows/RERCONVERGE/rer_trees"
include { RER_MATRIX }         from "${baseDir}/subworkflows/RERCONVERGE/rer_matrix"
include { RER_CONT }           from "${baseDir}/subworkflows/RERCONVERGE/rer_cont"
include { RER_BIN }            from "${baseDir}/subworkflows/RERCONVERGE/rer_bin"
include { RER_REPORT as RER_REPORT_CONT; RER_REPORT as RER_REPORT_BIN } from "${baseDir}/subworkflows/RERCONVERGE/rer_report.nf"

// Main workflow
workflow RER_MAIN {

    take:
        traitfile_input   // Channel<path> or null → falls back to params.my_traits
        // Kept for backwards-compatible signature; currently unused.
        pp_top_ch
        pp_bottom_ch

    main:

        // ── Resolve traitfile ────────────────────────────────────────────────
        def my_traitfile_ch = (traitfile_input ?: Channel.empty())
            .ifEmpty {
                assert params.my_traits : "RER_MAIN requires --my_traits"
                file(params.my_traits)
            }

        // ── Resolve gene trees (always full set) ─────────────────────────────
        def effective_gene_trees_ch = Channel.value(file(params.gene_trees))

        // ── Conditionally run RER steps ──────────────────────────────────────
        def trait_out_ch    = params.trait_out  ? Channel.value(file(params.trait_out))  : Channel.empty()
        def masterTrees_out = params.trees_out  ? Channel.value(file(params.trees_out))  : Channel.empty()
        def matrix_out_ch   = params.matrix_out ? Channel.value(file(params.matrix_out)) : Channel.empty()

        // trait_type_ch: string channel ("binary" or "continuous") resolved after
        // build_trait, or defaulting to params.rer_trait_mode when skipping build_trait.
        def trait_type_ch = Channel.empty()

        // Track RER_REPORT output for downstream consumers (e.g. SCORING).
        // Populated when 'continuous' or 'binary' tool is run; empty otherwise.
        def rer_report_out = null

        if (params.rer_tool) {
            def toolsToRun = params.rer_tool.split(',')

            if (toolsToRun.contains('build_trait')) {
                def trait_result = RER_TRAIT(my_traitfile_ch)
                trait_out_ch   = trait_result.polished
                trait_type_ch  = trait_result.trait_type.map { f -> f.text.trim() }
                trait_type_ch.view { t -> "[RER] Auto-detected trait type: ${t}" }
            }

            if (toolsToRun.contains('build_tree')) {
                def tax_id_ch = params.tax_id
                    ? Channel.value(file(params.tax_id))
                    : Channel.value(file('NO_FILE'))
                def trees_out = RER_TREES(my_traitfile_ch, effective_gene_trees_ch, tax_id_ch)
                masterTrees_out = trees_out[0]
            }

            if (toolsToRun.contains('build_matrix')) {
                matrix_out_ch = RER_MATRIX(trait_out_ch, masterTrees_out)
            }

            // ── Auto-route: continuous vs binary ─────────────────────────────
            // When rer_tool contains 'continuous' or 'binary', the resolved trait
            // type (from build_trait auto-detection or params.rer_trait_mode)
            // determines which correlation process runs.
            // Override auto-detection by setting params.rer_trait_mode = 'continuous'
            // or 'binary' explicitly.
            if (toolsToRun.contains('continuous') || toolsToRun.contains('binary')) {

                // Resolve effective type: explicit param overrides auto-detection.
                // If build_trait was skipped (trait_type_ch is empty) and no explicit
                // mode is set, default to 'continuous' for backwards compatibility.
                def effective_type_ch = (params.rer_trait_mode && params.rer_trait_mode != 'auto')
                    ? Channel.value(params.rer_trait_mode)
                    : trait_type_ch.ifEmpty { 'continuous' }

                // Combine trait file with its type string for branching
                def routed_trait = trait_out_ch
                    .combine(effective_type_ch)

                routed_trait.branch {
                    continuous: it[1] == 'continuous'
                    binary:     it[1] == 'binary'
                }.set { trait_branched }

                def cont_trait_ch = trait_branched.continuous.map { polished, _type -> polished }
                def bin_trait_ch  = trait_branched.binary.map     { polished, _type -> polished }

                def gmt_ch = params.rer_gmt_file
                    ? Channel.value(file(params.rer_gmt_file))
                    : Channel.value(file('NO_FILE'))

                // Both processes are declared; only the branch that received data
                // will produce output — the other stays idle (empty input channel).
                def rer_cont_result = RER_CONT(cont_trait_ch, masterTrees_out, matrix_out_ch)
                def rer_bin_result  = RER_BIN(bin_trait_ch,  masterTrees_out, matrix_out_ch)

                def cont_report = RER_REPORT_CONT(rer_cont_result.continuous_output, gmt_ch)
                def bin_report  = RER_REPORT_BIN(rer_bin_result.binary_output,      gmt_ch)

                // Merge summary TSV outputs (only the active branch will emit)
                rer_report_out = cont_report.summary_tsv.mix(bin_report.summary_tsv)
            }
        }

    emit:
        summary_tsv = rer_report_out ?: Channel.empty()
}
