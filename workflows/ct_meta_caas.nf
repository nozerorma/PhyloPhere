#!/usr/bin/env nextflow

/*
 * CT meta-CAAS workflow
 * Renders the CT pattern annotation report of a discovery.tab. The meta_caas tables are written by the observed
 * scoring (CAAS_CORE_OBSERVED, CAAS_OBSERVED).
 */

include { CAAS_META_CAAS_REPORT } from '../subworkflows/CT_META_CAAS/ctpp_meta_caas'

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
                global_background_genes = Channel.fromPath("${params.background_input}/*background_genes*", type: 'file')
                    .toList()
                    .map { files ->
                        if (files.size() != 1) {
                            error "CT_META_CAAS requires exactly one background_genes file in ${params.background_input}; found ${files.size()}: ${files}"
                        }
                        files[0]
                    }
                log.info "📂 Loading global background genes from directory: ${params.background_input}"
            } else {
                assert bg_path.isFile() : "CT_META_CAAS background_input must be a regular file or directory"
                global_background_genes = Channel.value(bg_path)
                log.info "📄 Loading global background genes file: ${params.background_input}"
            }
        } else {
            error "CT_META_CAAS requires CT background_genes output or --background_input"
        }

        // Header-only discovery is a valid result; the report renders its empty state.
        meta_caas_results = CAAS_META_CAAS_REPORT(
            discovery_file_ch,
            global_background_genes
        )

    emit:
        report = meta_caas_results.report
}
