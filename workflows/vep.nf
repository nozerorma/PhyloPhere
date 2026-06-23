#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: VEP characterization workflow
 *
 * Annotates CAAS positions with PrimateAI-3D pathogenicity scores.
 * Uses upstream MAP files directly.
 */

include { PRIMATEAI_MAP } from "${baseDir}/subworkflows/VEP/primateai.nf"

workflow VEP {
    take:
        caas_input

    main:
        // Output channels default to empty when VEP is enabled without any CAAS source.
        def primateai_out = Channel.empty()

        // Resolve CAAS input: integrated runs pass a channel; standalone runs use
        // --vep_caas_input.
        def caas_source = null
        if (caas_input) {
            caas_source = caas_input
        } else if (params.vep_caas_input) {
            def f = file(params.vep_caas_input)
            assert f.exists() : "VEP: CAAS input not found: ${params.vep_caas_input}"
            caas_source = Channel.value(f)
        }

        if (caas_source) {
            // .collect() + .map() converts to a value channel so the file can be
            // forwarded to downstream processes.
            def caas_ch = caas_source
                .collect()
                .filter { files -> files && files.size() > 0 }
                .map { files -> files[0] }

            // MAP files directory (upstream)
            assert params.vep_map_dir : "VEP requires --vep_map_dir (directory containing per-gene MAP TSV files)"
            def map_dir_ch = Channel.value(file(params.vep_map_dir))

            // Static reference data (defaults point to dat/ inside this subworkflow)
            def pai_db_ch  = Channel.value(file(params.vep_primateai_db))

            // ── PrimateAI-3D score mapping ───────────────────────────
            def pai_out = PRIMATEAI_MAP(caas_ch, map_dir_ch, pai_db_ch)

            primateai_out = pai_out.primateai_tsv
        } else {
            log.warn "VEP requested but no CAAS input was available from CT_POSTPROC and --vep_caas_input was not provided. Skipping VEP."
        }

    emit:
        primateai_tsv = primateai_out
}
