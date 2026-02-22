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
# File: ora_accumulation.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ORA_ACCUMULATION Workflow: Runs WebGestalt ORA (and optionally STRING) on gene lists
 *  produced by CT_ACCUMULATION_RANDOMIZE (accumulation/randomization/gene_lists/).
 *
 *  Inputs (dual-mode: integrated pipeline OR standalone params):
 *    - background_input_channel : cleaned_background_main.txt from CT_POSTPROC / CT_ACCUMULATION
 *    - gene_lists_input_channel : gene_lists channel emitted by CT_ACCUMULATION
 *
 *  Standalone fallbacks:
 *    - params.ora_background_input            : path to background .txt file or directory
 *    - params.accumulation_ora_gene_lists_input: path to directory containing gene list .txt files
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { ORA_ACCUMULATION_REPORT }    from "${baseDir}/subworkflows/ORA/ora_accumulation"
include { STRING_ACCUMULATION_REPORT } from "${baseDir}/subworkflows/ORA/string_accumulation"

workflow ORA_ACCUMULATION {
    take:
        background_input_channel   // cleaned background from CT_POSTPROC (or null for standalone)
        gene_lists_input_channel   // gene_lists from CT_ACCUMULATION (or null for standalone)

    main:
        // ── Background ────────────────────────────────────────────────────────
        def background_file_ch = (background_input_channel ?: Channel.empty())
        if (params.ora_background_input) {
            background_file_ch = background_file_ch.ifEmpty {
                def bg_path = file(params.ora_background_input)
                assert bg_path.exists() : "Error: ora_background_input not found: ${params.ora_background_input}"

                if (bg_path.isDirectory()) {
                    def preferred = file("${params.ora_background_input}/cleaned_background_main.txt")
                    if (preferred.exists()) {
                        log.info "📂 ORA_ACCUMULATION: using cleaned_background_main.txt from directory"
                        preferred
                    } else {
                        def txts = bg_path.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                        assert txts.size() > 0 : "Error: No .txt background files found in ${params.ora_background_input}"
                        log.info "📂 ORA_ACCUMULATION: using background file: ${txts[0]}"
                        txts[0]
                    }
                } else {
                    log.info "📄 ORA_ACCUMULATION: using background file: ${bg_path}"
                    bg_path
                }
            }
        }

        // ── Gene lists ────────────────────────────────────────────────────────
        // IMPORTANT:
        // - We intentionally avoid `if (gene_lists_input_channel) { ...collect()... }`
        //   because an upstream optional/empty channel can lead to a silent stall.
        // - Instead, we always build from the channel, then fallback with `.ifEmpty {}`.
        def gene_list_files_ch = (gene_lists_input_channel ?: Channel.empty())
            // CT_ACCUMULATION may emit gene_lists as either individual paths
            // or a single List<Path> item depending on glob expansion.
            .flatMap { item ->
                item instanceof Collection ? item : [item]
            }
            .filter { it.toString().endsWith('.txt') }
        if (params.accumulation_ora_gene_lists_input) {
            gene_list_files_ch = gene_list_files_ch.ifEmpty {
                def lists_dir = file(params.accumulation_ora_gene_lists_input)
                assert lists_dir.exists()    : "Error: accumulation_ora_gene_lists_input not found: ${params.accumulation_ora_gene_lists_input}"
                assert lists_dir.isDirectory(): "Error: accumulation_ora_gene_lists_input must be a directory"

                def txts = lists_dir.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                assert txts.size() > 0 : "Error: No .txt gene lists found in ${params.accumulation_ora_gene_lists_input}"

                log.info "📂 ORA_ACCUMULATION: using gene lists from directory: ${params.accumulation_ora_gene_lists_input} (${txts.size()} files)"
                txts
            }
        }

        gene_list_files_ch = gene_list_files_ch
            .collect()
            .filter { files -> files && files.flatten().findAll { it.toString().endsWith('.txt') }.size() > 0 }
            .map { files ->
                def resolved = files.flatten().findAll { it.toString().endsWith('.txt') }
                log.info "📥 ORA_ACCUMULATION: using ${resolved.size()} gene list(s)"
                resolved
            }

        ora_run = ORA_ACCUMULATION_REPORT(
            background_file_ch,
            gene_list_files_ch
        )

        // Run STRING enrichment in parallel when --string is enabled
        string_run = params.string
            ? STRING_ACCUMULATION_REPORT(background_file_ch, gene_list_files_ch)
            : null

    emit:
        report         = ora_run.report
        ora_results    = ora_run.ora_results
        ora_summary    = ora_run.ora_summary
        ora_plots      = ora_run.ora_plots
        string_report  = params.string ? string_run.report         : Channel.empty()
        string_summary = params.string ? string_run.string_summary : Channel.empty()
        string_plots   = params.string ? string_run.string_plots   : Channel.empty()
}
