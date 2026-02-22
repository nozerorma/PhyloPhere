#!/usr/bin/env nextflow

/*
 * CT disambiguation workflow
 * Runs the CT convergence-type disambiguation stage using metadata + ASR.
 */

include { CT_DISAMBIGUATION_RUN } from "${baseDir}/subworkflows/CT_DISAMBIGUATION/ct_disambiguation"

workflow CT_DISAMBIGUATION {
    take:
        meta_caas_in
        trait_file_in
        tree_file_in

    main:
        // Disambiguation must consume signification output from meta_caas/global_meta_caas.tsv
        def meta_from_upstream = (meta_caas_in ?: Channel.empty())
            .flatten()
            .filter { f ->
                def p = f.toString().toLowerCase()
                p.endsWith('global_meta_caas.tsv') || p.contains('meta_caas/global_meta_caas.tsv')
            }
            .collect()
            .map { files ->
                files ? files[0] : null
            }
            .filter { it != null }

        def meta_caas = meta_from_upstream.ifEmpty {
            if (params.ct_disambig_caas_metadata) {
                def f = file(params.ct_disambig_caas_metadata)
                assert f.exists() : "Error: ct_disambig_caas_metadata file not found: ${params.ct_disambig_caas_metadata}"
                log.info "📄 Loading standalone disambiguation metadata: ${params.ct_disambig_caas_metadata}"
                return f
            }
            error "CT disambiguation requires signification meta file or --ct_disambig_caas_metadata"
        }

        // Prefer upstream CT-resolved trait/tree for integrated runs; fall back to params for standalone runs
        def upstream_trait = (trait_file_in ?: Channel.empty())
        def upstream_tree = (tree_file_in ?: Channel.empty())

        def trait_file = upstream_trait.ifEmpty {
            def trait_file_param = params.caas_config
            if (!trait_file_param) {
                error "CT disambiguation requires a trait file from CT/contrast_selection or --caas_config"
            }
            file(trait_file_param)
        }

        def tree_file = upstream_tree.ifEmpty {
            def tree_file_param = params.tree
            if (!tree_file_param) {
                error "CT disambiguation requires a tree file from CT/contrast_selection or --tree"
            }
            file(tree_file_param)
        }

        CT_DISAMBIGUATION_RUN(meta_caas, trait_file, tree_file)

    emit:
        results_dir = CT_DISAMBIGUATION_RUN.out.results_dir
        master_csv = CT_DISAMBIGUATION_RUN.out.master_csv
}
