#!/usr/bin/env nextflow

/*
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: ora_excluded.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ORA Excluded Workflow: Runs WebGestalt ORA on excluded gene lists using the original
 *  (pre-cleanup) background from CT discovery.  Outputs to ${params.outdir}/ora_excluded
 *  so results are kept separate from the standard significant-gene ORA.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { ORA_EXCLUDED_REPORT } from "${baseDir}/subworkflows/ORA/ora_excluded_general"

workflow ORA_EXCLUDED {
    take:
        background_input_channel    // original (pre-cleanup) background gene file
        gene_lists_input_channel    // excluded gene list .txt files from CT_POSTPROC

    main:
        // Resolve background source: upstream channel OR standalone parameter.
        def background_file_ch
        if (background_input_channel) {
            background_file_ch = background_input_channel
        } else {
            assert params.ora_background_input : "ORA_EXCLUDED requires CT_POSTPROC background_ori output or --ora_background_input"

            def bg_path = file(params.ora_background_input)
            assert bg_path.exists() : "Error: ora_background_input not found: ${params.ora_background_input}"

            if (bg_path.isDirectory()) {
                def preferred = file("${params.ora_background_input}/background_genes.output")
                if (preferred.exists()) {
                    log.info "📂 ORA_EXCLUDED: using background file: ${preferred}"
                    background_file_ch = Channel.value(preferred)
                } else {
                    def txts = bg_path.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                    assert txts.size() > 0 : "Error: No .txt background files found in ${params.ora_background_input}"
                    background_file_ch = Channel.value(txts[0])
                }
            } else {
                background_file_ch = Channel.value(bg_path)
            }
        }

        // Resolve excluded gene-list source: upstream channel OR standalone parameter.
        def gene_list_files_ch
        if (gene_lists_input_channel) {
            gene_list_files_ch = gene_lists_input_channel
                .filter { it.toString().endsWith('.txt') }
                .collect()
                .map { files ->
                    log.info "📥 ORA_EXCLUDED: using ${files.size()} excluded gene list(s) from CT_POSTPROC"
                    files
                }
        } else {
            assert params.ora_gene_lists_input : "ORA_EXCLUDED requires CT_POSTPROC excluded gene-list outputs or --ora_gene_lists_input"

            def lists_dir = file(params.ora_gene_lists_input)
            assert lists_dir.exists()     : "Error: ora_gene_lists_input not found: ${params.ora_gene_lists_input}"
            assert lists_dir.isDirectory(): "Error: ora_gene_lists_input must be a directory"

            def txts = lists_dir.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
            assert txts.size() > 0 : "Error: No .txt gene lists found in ${params.ora_gene_lists_input}"

            log.info "📂 ORA_EXCLUDED: loading excluded gene lists from directory: ${params.ora_gene_lists_input} (${txts.size()} files)"
            gene_list_files_ch = Channel.value(txts)
        }

        ora_excl_run = ORA_EXCLUDED_REPORT(
            background_file_ch,
            gene_list_files_ch
        )

    emit:
        report      = ora_excl_run.report
        ora_results = ora_excl_run.ora_results
        ora_summary = ora_excl_run.ora_summary
        ora_plots   = ora_excl_run.ora_plots
}
