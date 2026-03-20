#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: CAAS Scoring Workflow
 *
 * Computes composite position-level and gene-level CAAS scores by
 * integrating outputs from CT_POSTPROC, PGLS, FADE, RERConverge,
 * and CT_ACCUMULATION.  Optionally runs ORA enrichment on top-scoring
 * gene lists.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 * File: workflows/scoring.nf
 */

include { SCORING_COMPUTE } from "${baseDir}/subworkflows/SCORING/scoring_compute.nf"
include { SCORING_REPORT }  from "${baseDir}/subworkflows/SCORING/scoring_report.nf"
include { ORA_GENERAL_REPORT as ORA_SCORING_POSITION } from "${baseDir}/subworkflows/ORA/ora_general.nf"
include { ORA_GENERAL_REPORT as ORA_SCORING_GENE     } from "${baseDir}/subworkflows/ORA/ora_general.nf"


workflow SCORING {

    take:
        postproc_ch              // Channel<path> or null — filtered_discovery.tsv
        pgls_ch                  // Channel<path> or null — site_pgls.tsv
        fade_summary_top_ch      // Channel<path> or null — fade_summary_top.tsv
        fade_summary_bottom_ch   // Channel<path> or null — fade_summary_bottom.tsv
        rer_summary_ch           // Channel<path> or null — rerconverge_summary_{trait}.tsv
        accum_ch                 // Channel<path> or null — collected accumulation CSVs
        background_ch            // Channel<path> or null — cleaned_background_main.txt

    main:
        assert params.traitname : "SCORING requires --traitname"

        // ── Resolve each input with .ifEmpty{} fallbacks ───────────────────

        def resolved_postproc = (postproc_ch ?: Channel.empty())
            .ifEmpty {
                assert params.scoring_postproc_input : \
                    "SCORING requires CT post-processing output or --scoring_postproc_input"
                def f = file(params.scoring_postproc_input)
                assert f.exists() : "SCORING: postproc input not found: ${params.scoring_postproc_input}"
                f
            }

        // Use unique sentinel names per input to avoid Nextflow staging collisions
        // (all file('NO_FILE') would try to stage as the same symlink name).
        def resolved_pgls = (pgls_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_pgls_input ?: 'NO_PGLS') }

        def resolved_fade_top = (fade_summary_top_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_summary_top ?: 'NO_FADE_TOP') }

        def resolved_fade_bottom = (fade_summary_bottom_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_summary_bottom ?: 'NO_FADE_BOTTOM') }

        def resolved_rer = (rer_summary_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_rer_input ?: 'NO_RER') }

        // Accumulation: if accum_ch is provided (value channel with collected files),
        // use it; otherwise check params.scoring_accum_dir for a directory of CSVs.
        def resolved_accum = (accum_ch ?: Channel.empty())
            .ifEmpty {
                def accum_dir = params.scoring_accum_dir ?: ''
                if (accum_dir && file(accum_dir).isDirectory()) {
                    // Collect all accumulation CSVs from the directory
                    def csvs = file(accum_dir).listFiles()
                        ?.findAll { it.name.endsWith('.csv') && it.name.startsWith('accumulation_') }
                    if (csvs && csvs.size() > 0) {
                        return csvs
                    }
                }
                file('NO_ACCUM')
            }

        // Background for ORA
        def resolved_background = (background_ch ?: Channel.empty())
            .ifEmpty {
                def bg = params.scoring_background_input ?: ''
                if (bg) {
                    def f = file(bg)
                    if (f.exists()) return f
                }
                file('NO_BACKGROUND')
            }

        // ── Run scoring computation ────────────────────────────────────────
        def compute_out = SCORING_COMPUTE(
            resolved_postproc,
            resolved_pgls,
            resolved_fade_top,
            resolved_fade_bottom,
            resolved_rer,
            resolved_accum
        )

        // ── Render report ──────────────────────────────────────────────────
        def report_out = SCORING_REPORT(
            compute_out.position_scores,
            compute_out.gene_scores,
            compute_out.gene_correlations
        )

        // ── ORA on scoring gene lists (optional) ───────────────────────────
        if (params.scoring_ora) {
            // Gate ORA behind the report: combine background with report output so
            // Nextflow only releases the background channel item after SCORING_REPORT
            // has finished, preventing ORA and SCORING_REPORT from racing over I/O.
            def bg_after_report = resolved_background
                .combine(report_out.report)
                .map { bg, _report -> bg }

            // Position-level ORA: gene lists from top-scoring positions
            ORA_SCORING_POSITION(
                bg_after_report,
                compute_out.position_gene_lists.collect()
            )

            // Gene-level ORA: gene lists from top-scoring genes
            ORA_SCORING_GENE(
                bg_after_report,
                compute_out.gene_gene_lists.collect()
            )
        }

    emit:
        position_scores = compute_out.position_scores
        gene_scores     = compute_out.gene_scores
        report          = report_out.report
}
