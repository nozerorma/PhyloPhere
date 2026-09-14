#!/usr/bin/env nextflow

/*
 * CT meta-CAAS workflow
 * Runs CT pattern annotation / meta_caas generation independently (upstream of disambiguation).
 */

include { CAAS_META_CAAS_REPORT } from "${baseDir}/subworkflows/CT_META_CAAS/ctpp_meta_caas"

workflow CT_META_CAAS {
    take:
        discovery_input_channel
        background_genes_channel

    main:
        // Discovery source
        def discovery_file_ch
        if (discovery_input_channel) {
            log.info "📥 Using discovery file from CT module for meta_caas"
            discovery_file_ch = discovery_input_channel
        } else {
            assert params.discovery_from : "CT_META_CAAS requires CT discovery output or --discovery_from"
            def discovery_file_obj = file(params.discovery_from)
            assert discovery_file_obj.exists() : "Error: discovery_from file not found: ${params.discovery_from}"
            assert discovery_file_obj.isFile() : "Error: discovery_from must be a file"
            discovery_file_ch = Channel.value(discovery_file_obj)
        }

        // Background genes source
        def global_background_genes
        if (background_genes_channel) {
            global_background_genes = background_genes_channel
            log.info "📥 Using global background genes from CT module for meta_caas"
        } else if (params.background_input) {
            def bg_path = file(params.background_input)
            assert bg_path.exists() : "Error: background_input file/directory not found: ${params.background_input}"

            if (bg_path.isDirectory()) {
                global_background_genes = Channel.fromPath("${params.background_input}/*background_genes*")
                log.info "📂 Loading global background genes from directory: ${params.background_input}"
            } else {
                global_background_genes = Channel.fromPath(params.background_input)
                log.info "📄 Loading global background genes file: ${params.background_input}"
            }
        } else {
            error "CT_META_CAAS requires CT background_genes output or --background_input"
        }

        // Guard: gracefully stop the pipeline when discovery has header only (no CAAS rows)
        def discovery_with_counts = discovery_file_ch
            .map { f ->
                def fh = (f instanceof java.nio.file.Path) ? f.toFile() : f
                def row_count = fh.readLines().size()
                tuple(f, row_count)
            }

        def discovery_file_nonempty = discovery_with_counts
            .filter { f, row_count ->
                def has_data = row_count > 1
                if (!has_data) {
                    exit 0, "No CAAS discoveries found (header-only discovery file: ${f}). Stopping pipeline gracefully."
                }
                return has_data
            }
            .map { f, row_count -> f }

        meta_caas_results = CAAS_META_CAAS_REPORT(
            discovery_file_nonempty,
            global_background_genes
        )

    emit:
        report = meta_caas_results.report
        meta_caas = meta_caas_results.meta_caas
        global_meta_caas = meta_caas_results.global_meta_caas
}
