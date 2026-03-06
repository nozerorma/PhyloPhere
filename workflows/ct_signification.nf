#!/usr/bin/env nextflow

/*
 * CT signification workflow
 * Runs CT signification independently (upstream of disambiguation).
 */

include { CAAS_SIGNIFICATION_REPORT } from "${baseDir}/subworkflows/CT_SIGNIFICATION/ctpp_signification"

workflow CT_SIGNIFICATION {
    take:
        discovery_input_channel
        background_genes_channel
        bootstrap_input_channel

    main:
        // Discovery source
        def discovery_file_ch
        if (discovery_input_channel) {
            log.info "📥 Using discovery file from CT module for signification"
            discovery_file_ch = discovery_input_channel
        } else {
            assert params.discovery_out : "CT signification requires CT discovery output or --discovery_out"
            def discovery_file_obj = file(params.discovery_out)
            assert discovery_file_obj.exists() : "Error: discovery_out file not found: ${params.discovery_out}"
            assert discovery_file_obj.isFile() : "Error: discovery_out must be a file"
            discovery_file_ch = Channel.value(discovery_file_obj)
        }

        // Background genes source
        def global_background_genes
        if (background_genes_channel) {
            global_background_genes = background_genes_channel
            log.info "📥 Using global background genes from CT module for signification"
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
            error "CT signification requires CT background_genes output or --background_input"
        }

        // Bootstrap source
        def bootstrap_files
        if (bootstrap_input_channel) {
            log.info "📥 Using bootstrap file from CT module"
            bootstrap_files = bootstrap_input_channel.map { file -> file }.collect()
        } else {
            assert params.bootstrap_input : "CT signification requires --bootstrap_input when not using CT pipeline"

            def bootstrap_path = file(params.bootstrap_input)
            assert bootstrap_path.exists() : "Error: bootstrap_input not found: ${params.bootstrap_input}"

            if (bootstrap_path.isDirectory()) {
                def boot_files = Channel.fromPath("${params.bootstrap_input}/*.{boot,tab}")
                bootstrap_files = boot_files.collect()
                def bootstrap_count = file(params.bootstrap_input).list().findAll {
                    it.endsWith('.boot') || it.endsWith('.tab')
                }.size()
                assert bootstrap_count > 0 : "Error: No .boot or .tab files found in directory ${params.bootstrap_input}"
                log.info "📂 Loading ${bootstrap_count} bootstrap files from directory: ${params.bootstrap_input}"
            } else {
                assert bootstrap_path.name.endsWith('.boot') || bootstrap_path.name.endsWith('.tab') :
                    "Error: bootstrap_input file must have .boot or .tab extension"
                bootstrap_files = Channel.fromPath(params.bootstrap_input).collect()
                log.info "📄 Loading single bootstrap file: ${params.bootstrap_input}"
            }
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

        signification_results = CAAS_SIGNIFICATION_REPORT(
            discovery_file_nonempty,
            global_background_genes,
            bootstrap_files
        )

    emit:
        signification_report = signification_results.report
        signification_gene_lists = signification_results.gene_lists
        signification_meta_caas = signification_results.meta_caas
        signification_global_meta = signification_results.global_meta_caas
}
