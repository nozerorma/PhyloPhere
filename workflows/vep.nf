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

include { EXTRACT_VEP_ARCHIVES } from "${baseDir}/subworkflows/VEP/vep_archives.nf"
include { PROT2COORD    } from "${baseDir}/subworkflows/VEP/prot2coord.nf"
include { PROT2AA       } from "${baseDir}/subworkflows/VEP/prot2aa.nf"
include { TRANSVAR_ANNO } from "${baseDir}/subworkflows/VEP/transvar.nf"
include { PRIMATEAI_MAP } from "${baseDir}/subworkflows/VEP/primateai.nf"

workflow VEP {
    take:
        caas_input

    main:
        def need_aa2nuc_stage  = !(params.vep_aa2nuc_input as String)
        def need_aa2prot_stage = !(params.vep_aa2prot_input as String)
        def need_alignment_dirs = need_aa2nuc_stage || need_aa2prot_stage

        if (need_alignment_dirs) {
            assert params.vep_cds_dir   : "VEP requires --vep_cds_dir (directory with per-gene CDS FASTA files) unless --vep_aa2nuc_input and --vep_aa2prot_input are provided"
            assert params.vep_track_dir : "VEP requires --vep_track_dir (directory with per-gene TRACK HTML files) unless --vep_aa2nuc_input and --vep_aa2prot_input are provided"
        }

        // Output channels default to empty when VEP is enabled without any CAAS source.
        def aa2nuc_out    = Channel.empty()
        def aa2prot_out   = Channel.empty()
        def transvar_out  = Channel.empty()
        def primateai_out = Channel.empty()

        // Resolve CAAS input: integrated runs pass a channel; standalone runs use
        // --vep_caas_input. main.nf passes null (not Channel.empty()) when there is
        // no upstream CAAS source so this fallback remains reachable.
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

            // Static reference data (defaults point to dat/ inside this subworkflow)
            def hs_cds_ch  = Channel.value(file(params.vep_hs_cds))
            def gff_ch     = Channel.value(file(params.vep_gff))
            def equiv_ch   = Channel.value(file(params.vep_gene_equiv))
            def tv_db_ch   = Channel.value(file(params.vep_transvar_db))
            def pai_db_ch  = Channel.value(file(params.vep_primateai_db))

            def vep_dirs = null
            if (need_alignment_dirs) {
                // Per-project directory inputs (project-specific, passed as params)
                def cds_dir_ch   = Channel.value(file(params.vep_cds_dir))
                def track_dir_ch = Channel.value(file(params.vep_track_dir))

                // ── Stage 0: extract per-gene files from archives (once, shared) ─
                vep_dirs = EXTRACT_VEP_ARCHIVES(caas_ch, cds_dir_ch, track_dir_ch)
            }

            // Optional precomputed inputs allow skipping Stage 1/2 while still
            // continuing with downstream TransVar and PrimateAI mapping.
            def aa2nuc_ch
            if (params.vep_aa2nuc_input) {
                def precompAa2nuc = file(params.vep_aa2nuc_input)
                assert precompAa2nuc.exists() : "VEP: precomputed aa2nuc input not found: ${params.vep_aa2nuc_input}"
                aa2nuc_ch = Channel.value(precompAa2nuc)
            } else {
                def coord_out = PROT2COORD(caas_ch, vep_dirs.cds_dir, vep_dirs.track_dir,
                                           hs_cds_ch, gff_ch, equiv_ch)
                aa2nuc_ch = coord_out.aa2nuc_csv
            }

            def aa2prot_ch
            if (params.vep_aa2prot_input) {
                def precompAa2prot = file(params.vep_aa2prot_input)
                assert precompAa2prot.exists() : "VEP: precomputed aa2prot input not found: ${params.vep_aa2prot_input}"
                aa2prot_ch = Channel.value(precompAa2prot)
            } else {
                def prot_out = PROT2AA(caas_ch, vep_dirs.cds_dir, vep_dirs.track_dir,
                                       hs_cds_ch, gff_ch, equiv_ch)
                aa2prot_ch = prot_out.aa2prot_csv
            }

            // ── Stage 3: TransVar annotation ──────────────────────────────────
            def tv_out = TRANSVAR_ANNO(aa2prot_ch, tv_db_ch)

            // ── Stage 4: PrimateAI-3D score mapping ───────────────────────────
            def pai_out = PRIMATEAI_MAP(tv_out.transvar_tsv,
                                        aa2prot_ch,
                                        caas_ch,
                                        pai_db_ch)

            aa2nuc_out    = aa2nuc_ch
            aa2prot_out   = aa2prot_ch
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
