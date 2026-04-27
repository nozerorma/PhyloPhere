#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: CAAS Scoring Workflow
 *
 * Computes position-level and gene-level CAAS scores by integrating outputs
 * from CT_POSTPROC, FADE, RERConverge, and CT_ACCUMULATION.
 * Runs once on the full postproc pool; directional characterisation is
 * performed post-scoring via the change_side column.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 * File: workflows/scoring.nf
 */

include { SCORING_COMPUTE }         from "${baseDir}/subworkflows/SCORING/scoring_compute.nf"
include { SCORING_REPORT }          from "${baseDir}/subworkflows/SCORING/scoring_report.nf"
include { ORA_GENERAL_REPORT as ORA_SCORING_POSITION   } from "${baseDir}/subworkflows/ORA/ora_general.nf"
include { ORA_GENERAL_REPORT as ORA_SCORING_GENE       } from "${baseDir}/subworkflows/ORA/ora_general.nf"
include { STRING_GENERAL_REPORT as STRING_SCORING_POSITION } from "${baseDir}/subworkflows/ORA/string_general.nf"
include { STRING_GENERAL_REPORT as STRING_SCORING_GENE     } from "${baseDir}/subworkflows/ORA/string_general.nf"


workflow SCORING {

    take:
        postproc_ch              // Channel<path> or null — filtered_discovery.tsv
        fade_summary_top_ch      // Channel<path> or null — fade_summary_top.tsv
        fade_summary_bottom_ch   // Channel<path> or null — fade_summary_bottom.tsv
        rer_summary_ch           // Channel<path> or null — rerconverge_summary_{trait}.tsv
        accum_ch                 // Channel<path> or null — collected accumulation CSVs (all directions)
        background_ch            // Channel<path> or null — cleaned_background_main.txt
        vep_transvar_ch          // Channel<path> or null — transvar_annotations.tsv (optional)
        vep_primateai_ch         // Channel<path> or null — primateai_scores.tsv     (optional)
        vep_aa2prot_ch           // Channel<path> or null — aa2prot_global.csv (optional)
        genomic_info_ch          // Channel<path> or null — gene genomic coords TSV  (optional)
        fade_site_top_ch         // Channel<path> or null — fade_site_bf_top.tsv     (optional)
        fade_site_bot_ch         // Channel<path> or null — fade_site_bf_bottom.tsv  (optional)

    main:
        assert params.traitname : "SCORING requires --traitname"

        // ── Resolve each input with .ifEmpty{} fallbacks ───────────────────

        def resolved_postproc
        if (postproc_ch) {
            resolved_postproc = postproc_ch
        } else if (params.scoring_postproc_input) {
            def f = file(params.scoring_postproc_input)
            assert f.exists() : "SCORING: postproc input not found: ${params.scoring_postproc_input}"
            resolved_postproc = Channel.value(f)
        } else {
            resolved_postproc = Channel.empty()
        }

        def resolved_fade_top = (fade_summary_top_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_summary_top ?: 'NO_FADE_TOP') }

        def resolved_fade_bottom = (fade_summary_bottom_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_summary_bottom ?: 'NO_FADE_BOTTOM') }

        def resolved_rer = (rer_summary_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_rer_input ?: 'NO_RER') }

        // Accumulation: collect all CSVs (top + bottom + all) into a single value channel.
        // The R script selects accumulation_all_* files via pattern matching.
        def accum_source = (accum_ch ?: Channel.empty())
        if (params.scoring_accum_dir) {
            def accum_dir = params.scoring_accum_dir ?: ''
            accum_source = accum_source.ifEmpty {
                if (accum_dir && file(accum_dir).isDirectory()) {
                    def csvs = file(accum_dir).listFiles()
                        ?.findAll { it.name.endsWith('.csv') && it.name.startsWith('accumulation_') }
                    if (csvs && csvs.size() > 0) return csvs
                }
                [file('NO_ACCUM')]
            }
        }
        def accum_all_ch = accum_source
            .flatten()
            .collect()
            .ifEmpty { [file('NO_ACCUM')] }

        // VEP optional inputs
        def resolved_vep_transvar = (vep_transvar_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_transvar ?: 'NO_VEP_TRANSVAR') }

        def resolved_vep_primateai = (vep_primateai_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_primateai ?: 'NO_VEP_PRIMATEAI') }

        def resolved_vep_aa2prot = (vep_aa2prot_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_aa2prot ?: 'NO_VEP_AA2PROT') }

        def resolved_genomic_info = (genomic_info_ch ?: Channel.empty())
            .ifEmpty {
                def gi = params.gene_ensembl_file ?: ''
                if (gi && file(gi).exists()) file(gi) else file('NO_GENOMIC_INFO')
            }

        def resolved_fade_site_top_ch = (fade_site_top_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_site_top ?: 'NO_FADE_SITE_TOP.txt') }

        def resolved_fade_site_bot_ch = (fade_site_bot_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_site_bottom ?: 'NO_FADE_SITE_BOT.txt') }

        def resolved_background = (background_ch ?: Channel.empty())
            .ifEmpty {
                def bg = params.scoring_background_input ?: ''
                if (bg) {
                    def f = file(bg)
                    if (f.exists()) return f
                }
                file('NO_BACKGROUND')
            }

        // ── Run scoring — single pass on full postproc pool ────────────────
        def compute_out = SCORING_COMPUTE(
            resolved_postproc,
            resolved_fade_top,
            resolved_fade_bottom,
            resolved_rer,
            accum_all_ch
        )

        // ── Render report ──────────────────────────────────────────────────
        // Helper: for optional compute outputs, emit a sentinel when absent.
        def _opt = { ch, sentinel ->
            ch.ifEmpty { file(sentinel) }
        }

        def report_out = SCORING_REPORT(
            compute_out.position_scores,
            compute_out.gene_scores,
            compute_out.gene_correlations,
            _opt(compute_out.stress_summary,         'NO_SCORING_STRESS_SUMMARY'),
            _opt(compute_out.stress_correlations,    'NO_SCORING_STRESS_CORR'),
            _opt(compute_out.stress_rank_agreement,  'NO_SCORING_STRESS_RANK'),
            _opt(compute_out.stress_top_overlap,     'NO_SCORING_STRESS_OVERLAP'),
            _opt(compute_out.stress_variants,        'NO_SCORING_STRESS_VARIANTS'),
            _opt(compute_out.stress_latent_loadings, 'NO_SCORING_STRESS_LOADINGS'),
            resolved_fade_site_top_ch,
            resolved_fade_site_bot_ch,
            resolved_vep_transvar,
            resolved_vep_primateai,
            resolved_vep_aa2prot,
            resolved_genomic_info
        )

        // ── ORA on scoring gene lists (optional) ───────────────────────────
        if (params.scoring_ora) {
            def bg_after_report = resolved_background
                .combine(report_out.report.collect())
                .map { it[0] }

            ORA_SCORING_POSITION(
                bg_after_report,
                compute_out.position_gene_lists.collect(),
                'ORA_scoring_position'
            )

            ORA_SCORING_GENE(
                bg_after_report,
                compute_out.gene_gene_lists.collect(),
                'ORA_scoring_gene'
            )

            if (params.string) {
                STRING_SCORING_POSITION(
                    bg_after_report,
                    compute_out.position_gene_lists.collect()
                )
                STRING_SCORING_GENE(
                    bg_after_report,
                    compute_out.gene_gene_lists.collect()
                )
            }
        }

    emit:
        position_scores = compute_out.position_scores
        gene_scores     = compute_out.gene_scores
        report          = report_out.report
}
