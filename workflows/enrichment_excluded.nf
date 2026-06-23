#!/usr/bin/env nextflow

/*
#  Enrichment Excluded Workflow: Runs clusterProfiler enrichment on excluded gene lists
#  using the original (pre-cleanup) background from CT discovery.
*/

include { EXCLUDED_ENRICHMENT_REPORT } from "${baseDir}/subworkflows/ENRICHMENT/excluded_enrichment"

workflow ENRICHMENT_EXCLUDED {
    take:
        background_input_channel    // original (pre-cleanup) background gene file
        gene_lists_input_channel    // excluded gene list .txt files from CT_POSTPROC

    main:
        // Resolve background source: upstream channel OR standalone parameter.
        def background_file_ch
        if (background_input_channel) {
            background_file_ch = background_input_channel
        } else {
            assert params.enrichment_background_input : "ENRICHMENT_EXCLUDED requires CT_POSTPROC background_ori output or --enrichment_background_input"

            def bg_path = file(params.enrichment_background_input)
            assert bg_path.exists() : "Error: enrichment_background_input not found: ${params.enrichment_background_input}"

            if (bg_path.isDirectory()) {
                def preferred = file("${params.enrichment_background_input}/background_genes.output")
                if (preferred.exists()) {
                    log.info "📂 ENRICHMENT_EXCLUDED: using background file: ${preferred}"
                    background_file_ch = Channel.value(preferred)
                } else {
                    def txts = bg_path.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
                    assert txts.size() > 0 : "Error: No .txt background files found in ${params.enrichment_background_input}"
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
                    log.info "📥 ENRICHMENT_EXCLUDED: using ${files.size()} excluded gene list(s) from CT_POSTPROC"
                    files
                }
        } else {
            assert params.enrichment_gene_lists_input : "ENRICHMENT_EXCLUDED requires CT_POSTPROC excluded gene-list outputs or --enrichment_gene_lists_input"

            def lists_dir = file(params.enrichment_gene_lists_input)
            assert lists_dir.exists()     : "Error: enrichment_gene_lists_input not found: ${params.enrichment_gene_lists_input}"
            assert lists_dir.isDirectory(): "Error: enrichment_gene_lists_input must be a directory"

            def txts = lists_dir.listFiles()?.findAll { it.isFile() && it.name.endsWith('.txt') } ?: []
            assert txts.size() > 0 : "Error: No .txt gene lists found in ${params.enrichment_gene_lists_input}"

            log.info "📂 ENRICHMENT_EXCLUDED: loading excluded gene lists from directory: ${params.enrichment_gene_lists_input} (${txts.size()} files)"
            gene_list_files_ch = Channel.value(txts)
        }

        ora_excl_run = EXCLUDED_ENRICHMENT_REPORT(
            background_file_ch,
            gene_list_files_ch
        )

    emit:
        report      = ora_excl_run.report
        ora_results = ora_excl_run.ora_results
        ora_summary = ora_excl_run.ora_summary
        ora_plots   = ora_excl_run.ora_plots
}
