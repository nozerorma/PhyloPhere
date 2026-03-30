#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: VEP characterization workflow
 *
 * Annotates CAAS positions with genomic coordinates, TransVar variant
 * annotations, and PrimateAI-3D pathogenicity scores.
 *
 * Stage order:
 *   prot2coord ──┐
 *                ├── (parallel) ──→ transvar ──→ primateai
 *   prot2aa    ──┘
 *
 * prot2coord and prot2aa run in parallel (identical inputs, independent outputs).
 * prot2aa feeds transvar, which feeds primateai.
 * prot2coord output (aa2nuc_global.csv) is a standalone characterization artifact.
 */

include { PROT2COORD    } from "${baseDir}/subworkflows/VEP/prot2coord.nf"
include { PROT2AA       } from "${baseDir}/subworkflows/VEP/prot2aa.nf"
include { TRANSVAR_ANNO } from "${baseDir}/subworkflows/VEP/transvar.nf"
include { PRIMATEAI_MAP } from "${baseDir}/subworkflows/VEP/primateai.nf"

workflow VEP {
    take:
        caas_input

    main:
        assert params.vep_cds_dir   : "VEP requires --vep_cds_dir (directory with per-gene CDS FASTA files)"
        assert params.vep_track_dir : "VEP requires --vep_track_dir (directory with per-gene TRACK HTML files)"

        // Output channels default to empty when VEP is enabled without any CAAS source.
        def aa2nuc_out    = Channel.empty()
        def aa2prot_out   = Channel.empty()
        def transvar_out  = Channel.empty()
        def primateai_out = Channel.empty()

        // Resolve CAAS input: integrated runs pass a channel; standalone runs use
        // --vep_caas_input. Do not assert on empty upstream channels here; if a
        // previous workflow aborts or produces no CAAS table, VEP should stay idle
        // instead of throwing a second, misleading error.
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
            // forwarded to multiple downstream processes.
            def caas_ch = caas_source
                .collect()
                .filter { files -> files && files.size() > 0 }
                .map { files -> files[0] }

            // Per-project directory inputs (project-specific, passed as params)
            def cds_dir_ch   = Channel.value(file(params.vep_cds_dir))
            def track_dir_ch = Channel.value(file(params.vep_track_dir))

            // Static reference data (defaults point to dat/ inside this subworkflow)
            def hs_cds_ch  = Channel.value(file(params.vep_hs_cds))
            def gff_ch     = Channel.value(file(params.vep_gff))
            def equiv_ch   = Channel.value(file(params.vep_gene_equiv))
            def tv_db_ch   = Channel.value(file(params.vep_transvar_db))
            def pai_db_ch  = Channel.value(file(params.vep_primateai_db))

            // ── Stage 1+2 (parallel) ──────────────────────────────────────────
            def coord_out = PROT2COORD(caas_ch, cds_dir_ch, track_dir_ch,
                                       hs_cds_ch, gff_ch, equiv_ch)
            def prot_out  = PROT2AA   (caas_ch, cds_dir_ch, track_dir_ch,
                                       hs_cds_ch, gff_ch, equiv_ch)

            // ── Stage 3: TransVar annotation ──────────────────────────────────
            def tv_out = TRANSVAR_ANNO(prot_out.aa2prot_csv, tv_db_ch)

            // ── Stage 4: PrimateAI-3D score mapping ───────────────────────────
            def pai_out = PRIMATEAI_MAP(tv_out.transvar_tsv,
                                        prot_out.aa2prot_csv,
                                        caas_ch,
                                        pai_db_ch)

            aa2nuc_out    = coord_out.aa2nuc_csv
            aa2prot_out   = prot_out.aa2prot_csv
            transvar_out  = tv_out.transvar_tsv
            primateai_out = pai_out.primateai_tsv
        } else {
            log.warn "VEP requested but no CAAS input was available from CT_POSTPROC and --vep_caas_input was not provided. Skipping VEP."
        }

    emit:
        aa2nuc_csv    = aa2nuc_out
        aa2prot_csv   = aa2prot_out
        transvar_tsv  = transvar_out
        primateai_tsv = primateai_out
}
