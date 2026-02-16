#!/usr/bin/env nextflow

/*
 * CT disambiguation workflow
 * Runs the CAAS convergence-type disambiguation stage using metadata + ASR.
 */

include { CT_DISAMBIGUATION_RUN } from "${baseDir}/subworkflows/CT_DISAMBIGUATION/ct_disambiguation"

workflow CT_DISAMBIGUATION {
    take:
        meta_caas_in

    main:
        def meta_from_params = Channel.value(file(params.ct_disambig_caas_metadata))

        // Upstream metadata can come as a file list (e.g. CAAS_SIGNIFICATION emit)
        // Keep only tabular files likely to contain meta-CAAS and take the first one.
        def meta_from_upstream = (meta_caas_in ?: Channel.empty())
            .flatten()
            .filter { f ->
                def p = f.toString().toLowerCase()
                p.endsWith('.tsv') || p.endsWith('.tab') || p.endsWith('.csv') || p.endsWith('.txt')
            }
            .first()

        def meta_caas = meta_from_upstream.ifEmpty { meta_from_params }

        def trait_file = Channel.value(file(params.ct_disambig_trait_file))
        def tree_file = Channel.value(file(params.ct_disambig_tree))

        CT_DISAMBIGUATION_RUN(meta_caas, trait_file, tree_file)

    emit:
        results_dir = CT_DISAMBIGUATION_RUN.out.results_dir
        master_csv = CT_DISAMBIGUATION_RUN.out.master_csv
}
