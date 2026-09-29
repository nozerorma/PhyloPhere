#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: VEP characterization workflow
 *
 * Annotates CAAS positions with PrimateAI-3D pathogenicity scores.
 * Uses upstream MAP files directly.
 */

include { PRIMATEAI_MAP } from '../subworkflows/VEP/primateai.nf'
include { COSMIC_MAP }     from '../subworkflows/VEP/cosmic.nf'
include { ENSEMBL_VEP_ANNOTATE } from '../subworkflows/VEP/ensembl_vep.nf'

workflow VEP {
    take:
        position_scores_input   // SCORING position_scores.tsv

    main:
        // Output channels default to empty when VEP is enabled without databases.
        def primateai_out = Channel.empty()
        def cosmic_out = Channel.empty()
        def ensembl_vep_out = Channel.empty()

        // Resolve position_scores channel
        def ps_ch = null
        if (position_scores_input) {
            ps_ch = position_scores_input
                .collect()
                .filter { files -> files && files.size() > 0 }
                .map { files -> files[0] }
        }

        if (ps_ch) {
            // MAP files directory (upstream)
            assert params.vep_map_dir : "VEP requires --vep_map_dir (directory containing per-gene MAP TSV files)"
            def map_dir_ch = Channel.value(file(params.vep_map_dir))

            // ── PrimateAI-3D score mapping (conditional on database existence) ──
            def pai_db_file = params.vep_primateai_db ? file(params.vep_primateai_db) : file('NO_FILE')
            if (pai_db_file.name != 'NO_FILE' && pai_db_file.exists() && pai_db_file.size() > 0) {
                def pai_db_ch = Channel.value(pai_db_file)
                def pai_out = PRIMATEAI_MAP(ps_ch, map_dir_ch, pai_db_ch)
                primateai_out = pai_out.primateai_tsv
            } else {
                log.info "ℹ VEP: PrimateAI-3D database not provided/empty — skipping PrimateAI-3D pathogenicity mapping."
            }

            // ── COSMIC mapping (conditional on database existence) ──
            def cosmic_db_file = params.cosmic_db ? file(params.cosmic_db) : file('NO_FILE')
            if (cosmic_db_file.name != 'NO_FILE' && cosmic_db_file.exists() && cosmic_db_file.size() > 0) {
                def cosmic_db_ch = Channel.value(cosmic_db_file)
                COSMIC_MAP(ps_ch, map_dir_ch, cosmic_db_ch)
                cosmic_out = COSMIC_MAP.out.cosmic_tsv
            } else {
                log.info "ℹ VEP: COSMIC database not provided/empty — skipping COSMIC somatic mutation mapping."
            }

            // ── Ensembl VEP consequence annotation (independent of the DBs above) ──
            if (params.vep_ensembl) {
                assert params.gene_ensembl_file : "VEP: --vep_ensembl requires --gene_ensembl_file (for human_protein_id -> HGVS lookup)."
                // Not user-required: an empty --vep_cache_dir resolves to a
                // persistent, species/assembly-scoped default that
                // ENSEMBL_VEP_ANNOTATE populates itself on first use via
                // vep_install (see subworkflows/VEP/ensembl_vep.nf).
                def species = params.vep_species ?: 'homo_sapiens'
                def assembly = params.vep_assembly ?: 'GRCh38'
                def resolved_cache_dir = params.vep_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/vep/${species}_${assembly}"
                def cache_dir_ch = Channel.value(resolved_cache_dir)
                def ensembl_file_ch = Channel.value(file(params.gene_ensembl_file))
                def ensembl_out = ENSEMBL_VEP_ANNOTATE(ps_ch, map_dir_ch, ensembl_file_ch, cache_dir_ch)
                ensembl_vep_out = ensembl_out.ensembl_vep_tsv
            }
        } else {
            log.warn "VEP requested but position_scores input was empty. Skipping VEP."
        }

    emit:
        primateai_tsv = primateai_out
        cosmic_tsv = cosmic_out
        ensembl_vep_tsv = ensembl_vep_out
}
