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
# File: ora.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ORA Workflow: Runs WebGestalt ORA from CT post-processing outputs (integrated mode)
 *  or user-provided background/list inputs (standalone mode).
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { ORA_GENERAL_REPORT } from "${baseDir}/subworkflows/ORA/ora_general"
include { STRING_GENERAL_REPORT } from "${baseDir}/subworkflows/ORA/string_general"

workflow ORA {
    take:
        background_input_channel
        gene_lists_input_channel

    main:
        // Resolve background source: upstream channel OR standalone parameter.
        // Use an upfront if(channel) guard (Groovy-time) instead of .ifEmpty{}
        // (reactive-time) so the fallback is not triggered prematurely when
        // background_input_channel is a value channel that has not yet emitted
        // at DAG-construction time (same pattern used by CT_POSTPROC, FADE, etc.).
        def background_file_ch
        if (background_input_channel) {
            background_file_ch = background_input_channel
        } else {
            assert params.ora_background_input : "ORA requires CT_POSTPROC cleaned background output or --ora_background_input"

            def bg_path = file(params.ora_background_input)
            assert bg_path.exists() : "Error: ora_background_input not found: ${params.ora_background_input}"

            if (bg_path.isDirectory()) {
                def preferred = file("${params.ora_background_input}/cleaned_background_main.txt")
                if (preferred.exists()) {
                    log.info "📂 Using preferred cleaned background file: ${preferred}"
                    background_file_ch = Channel.value(preferred)
                } else {
                    def txts = bg_path.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                    assert txts.size() > 0 : "Error: No .txt background files found in ${params.ora_background_input}"
                    log.info "📂 Using first background .txt file from directory: ${txts[0]}"
                    background_file_ch = Channel.value(txts[0])
                }
            } else {
                log.info "📄 Using ORA background file: ${bg_path}"
                background_file_ch = Channel.value(bg_path)
            }
        }

        // Resolve gene-list source: upstream channel OR standalone parameter.
        // IMPORTANT: Channel.empty().collect() never emits in Nextflow DSL2 (it silently
        // completes without producing an item), so .map{} fallbacks after .collect() are
        // unreliable for empty channels. Use an upfront if(channel) guard instead,
        // matching the same pattern used by CT_POSTPROC for its disambiguation_input.
        def gene_list_files_ch
        if (gene_lists_input_channel) {
            gene_list_files_ch = gene_lists_input_channel
                .filter { it.toString().endsWith('.txt') }
                .collect()
                .map { files ->
                    log.info "📥 Using ${files.size()} gene list(s) from upstream CT_POSTPROC output"
                    files
                }
        } else {
            assert params.ora_gene_lists_input : "ORA requires CT_POSTPROC gene-list outputs or --ora_gene_lists_input"

            def lists_dir = file(params.ora_gene_lists_input)
            assert lists_dir.exists() : "Error: ora_gene_lists_input not found: ${params.ora_gene_lists_input}"
            assert lists_dir.isDirectory() : "Error: ora_gene_lists_input must be a directory"

            def txts = lists_dir.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
            assert txts.size() > 0 : "Error: No .txt gene lists found in ${params.ora_gene_lists_input}"

            log.info "📂 Using ORA gene lists from directory: ${params.ora_gene_lists_input} (${txts.size()} files)"
            gene_list_files_ch = Channel.value(txts)
        }

        ora_run = ORA_GENERAL_REPORT(
            background_file_ch,
            gene_list_files_ch
        )

        // Run STRING enrichment in parallel when --string is enabled
        string_run = params.string
            ? STRING_GENERAL_REPORT(background_file_ch, gene_list_files_ch)
            : null

    emit:
        report      = ora_run.report
        ora_results = ora_run.ora_results
        ora_summary = ora_run.ora_summary
        ora_plots   = ora_run.ora_plots
        string_report  = params.string ? string_run.report        : Channel.empty()
        string_summary = params.string ? string_run.string_summary : Channel.empty()
        string_plots   = params.string ? string_run.string_plots   : Channel.empty()
}


