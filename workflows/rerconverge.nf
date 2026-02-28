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
include { COLLECT_GENE_SETS }  from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { FILTER_GENE_TREES }  from "${baseDir}/subworkflows/RERCONVERGE/rer_filter_trees.nf"

// Main workflow
workflow RER_MAIN {

    take:
        traitfile_input   // Channel<path> or null → falls back to params.my_traits
        // Optional upstream channels for inline gene-set piping.
        // Pass Channel.empty() when running standalone (params-based).
        acc_top_ch        // accumulation CSV for TOP
        acc_bottom_ch     // accumulation CSV for BOTTOM
        pp_top_ch         // postproc TXT for TOP
        pp_bottom_ch      // postproc TXT for BOTTOM

    main:

        // ── Resolve traitfile ────────────────────────────────────────────────
        def my_traitfile_ch = (traitfile_input ?: Channel.empty())
            .ifEmpty {
                assert params.my_traits : "RER_MAIN requires --my_traits"
                file(params.my_traits)
            }

        // ── Resolve gene trees (optionally pre-filtered to gene set) ─────────
        def effective_gene_trees_ch

        if (params.rer_gene_set_mode == 'gene_set') {
            // Prefer upstream piped channels; fall back to --rer_* path params.
            def resolved_acc_top = (acc_top_ch ?: Channel.empty())
                .ifEmpty { file(params.rer_accumulation_top    ?: 'NO_FILE') }

            def resolved_acc_bottom = (acc_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.rer_accumulation_bottom ?: 'NO_FILE') }

            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.rer_postproc_top        ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.rer_postproc_bottom     ?: 'NO_FILE') }

            def gene_sets = COLLECT_GENE_SETS(
                resolved_acc_top,
                resolved_acc_bottom,
                resolved_pp_top,
                resolved_pp_bottom
            )

            def filtered = FILTER_GENE_TREES(
                Channel.value(file(params.gene_trees)),
                gene_sets.gene_set_top,
                gene_sets.gene_set_bottom
            )
            effective_gene_trees_ch = filtered.filtered_trees

        } else {
            // all mode: use the full gene trees file unchanged
            effective_gene_trees_ch = Channel.value(file(params.gene_trees))
        }

        // ── Conditionally run RER steps ──────────────────────────────────────
        def trait_out      = params.trait_out     ? Channel.value(file(params.trait_out))     : Channel.empty()
        def masterTrees_out = params.trees_out    ? Channel.value(file(params.trees_out))     : Channel.empty()
        def matrix_out_ch  = params.matrix_out    ? Channel.value(file(params.matrix_out))    : Channel.empty()

        if (params.rer_tool) {
            def toolsToRun = params.rer_tool.split(',')

            if (toolsToRun.contains('build_trait')) {
                trait_out = RER_TRAIT(my_traitfile_ch)
            }

            if (toolsToRun.contains('build_tree')) {
                trees_out = RER_TREES(my_traitfile_ch, effective_gene_trees_ch)
                masterTrees_out = trees_out[0]
            }

            if (toolsToRun.contains('build_matrix')) {
                matrix_out_ch = RER_MATRIX(trait_out, masterTrees_out)
            }

            if (toolsToRun.contains('continuous')) {
                RER_CONT(trait_out, masterTrees_out, matrix_out_ch)
            }
        }
}
