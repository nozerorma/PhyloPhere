#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: VEP characterization workflow
 *
 * Annotates CAAS positions with PrimateAI-3D pathogenicity scores.
 * Uses upstream MAP files directly.
 */

include { PRIMATEAI_MAP } from "${baseDir}/subworkflows/VEP/primateai.nf"
include { COSMIC_MAP }     from "${baseDir}/subworkflows/VEP/cosmic.nf"

workflow VEP {
    take:
        caas_input

    main:
        // Output channels default to empty when VEP is enabled without any CAAS source.
        def primateai_out = Channel.empty()
        def cosmic_out = Channel.empty()

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

            // ── PrimateAI-3D score mapping (conditional on database existence) ──
            def pai_db_file = params.vep_primateai_db ? file(params.vep_primateai_db) : file('NO_FILE')
            if (pai_db_file.name != 'NO_FILE' && pai_db_file.exists() && pai_db_file.size() > 0) {
                def pai_db_ch = Channel.value(pai_db_file)
                def pai_out = PRIMATEAI_MAP(caas_ch, map_dir_ch, pai_db_ch)
                primateai_out = pai_out.primateai_tsv
            } else {
                log.info "ℹ VEP: PrimateAI-3D database not provided/empty — skipping PrimateAI-3D pathogenicity mapping."
            }

            // ── COSMIC mapping (conditional on database existence) ──
            def cosmic_db_file = params.cosmic_db ? file(params.cosmic_db) : file('NO_FILE')
            if (cosmic_db_file.name != 'NO_FILE' && cosmic_db_file.exists() && cosmic_db_file.size() > 0) {
                def cosmic_db_ch = Channel.value(cosmic_db_file)
                COSMIC_MAP(caas_ch, map_dir_ch, cosmic_db_ch)
                cosmic_out = COSMIC_MAP.out.cosmic_tsv
            } else {
                log.info "ℹ VEP: COSMIC database not provided/empty — skipping COSMIC somatic mutation mapping."
            }
        } else {
            log.warn "VEP requested but no CAAS input was available from CT_POSTPROC and --vep_caas_input was not provided. Skipping VEP."
        }

    emit:
        primateai_tsv = primateai_out
        cosmic_tsv = cosmic_out
}
