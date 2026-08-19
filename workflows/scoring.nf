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

include { SCORING_COMPUTE }                                                           from "${baseDir}/subworkflows/SCORING/scoring_compute.nf"
include { SCORING_REPORT }                                                            from "${baseDir}/subworkflows/SCORING/scoring_report.nf"
include { CAAS_PERMS_REBUILD }                                                        from "${baseDir}/subworkflows/CT/caas_permulation.nf"


workflow SCORING {

    take:
        postproc_ch              // Channel<path> or null — filtered_discovery.tsv
        fade_summary_top_ch      // Channel<path> or null — fade_summary_top.tsv
        fade_summary_bottom_ch   // Channel<path> or null — fade_summary_bottom.tsv
        rer_summary_ch           // Channel<path> or null — rerconverge_summary_{trait}.tsv
        accum_ch                 // Channel<path> or null — collected accumulation CSVs (all directions)
        genomic_info_ch          // Channel<path> or null — gene genomic coords TSV  (optional)
        fade_site_top_ch         // Channel<path> or null — fade_site_bf_top.tsv     (optional)
        fade_site_bot_ch         // Channel<path> or null — fade_site_bf_bottom.tsv  (optional)
        cleaned_background_ch    // Channel<path> or null — cleaned_background_main.txt (FCS universe)
        rer_perms_ch             // Channel<path> or null — RER permulation RDS (corStat) for RER FCS p.perm
        caas_perms_ch            // Channel<path> or null — CAAS permulation RDS (asr + caas null) for FCS p.perm + report
        caas_pos_pval_ch         // Channel<path> or null — LOO null_pvalue_boot per (gene,position,scheme)
        caas_pos_sample_ch       // Channel<path> or null — cycle-stratified per-scheme sample for report distribution plots
        caas_pos_quantiles_ch    // Channel<path> or null — per (cycle,scheme) null distribution shape

    main:
        assert params.traitname : "SCORING requires --traitname"

        // ── Resolve each input with .ifEmpty{} fallbacks ───────────────────

        def resolved_postproc
        if (postproc_ch) {
            // .first() → broadcastable value channel: filtered_discovery is consumed
            // by BOTH SCORING_COMPUTE and the report's functional-group overlay tab.
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
            .ifEmpty { file(params.scoring_rer_input ?: (params.rer_continuous_file ?: 'NO_RER')) }

        // Accumulation: collect all CSVs (top + bottom + all) into a single value channel.
        // The R script selects accumulation_all_* files via pattern matching.
        def accum_source = (accum_ch ?: Channel.empty())
        if (params.scoring_accum_dir) {
            def accum_dir = params.scoring_accum_dir ?: ''
            accum_source = accum_source.ifEmpty {
                // Recursive: live output is accumulation/{top,bottom,all}/randomization/*.csv,
                // three levels below accum_dir, not flat inside it — a plain listFiles() here
                // would find zero matches against the real directory-layout.
                if (accum_dir && file(accum_dir).isDirectory()) {
                    def csvs = []
                    file(accum_dir).eachFileRecurse(groovy.io.FileType.FILES) { f ->
                        if (f.name.endsWith('.csv') && f.name.startsWith('accumulation_')) csvs << f
                    }
                    if (csvs.size() > 0) return csvs
                }
                [file('NO_ACCUM')]
            }
        }
        def accum_all_ch = accum_source
            .flatten()
            .collect()
            .ifEmpty { [file('NO_ACCUM')] }

        def resolved_genomic_info = (genomic_info_ch ?: Channel.empty())
            .ifEmpty {
                def gi = params.gene_ensembl_file ?: ''
                if (gi && file(gi).exists()) file(gi) else file('NO_GENOMIC_INFO')
            }

        // .collect() turns these into value channels so they can be safely consumed
        // by both SCORING_COMPUTE and SCORING_REPORT below — without it, a regular
        // (queue) channel is drained by its first consumer, leaving the report with
        // nothing and silently dropping the FADE site data from the rendered report.
        // .map { it[0] } unwraps the single-element list .collect() produces back to
        // a scalar path, since both processes call `.name` directly on this input.
        def resolved_fade_site_top_ch = (fade_site_top_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_site_top ?: 'NO_FADE_SITE_TOP.txt') }
            .collect()
            .map { it && it.size() > 0 ? it[0] : file('NO_FADE_SITE_TOP.txt') }

        def resolved_fade_site_bot_ch = (fade_site_bot_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_fade_site_bottom ?: 'NO_FADE_SITE_BOT.txt') }
            .collect()
            .map { it && it.size() > 0 ? it[0] : file('NO_FADE_SITE_BOT.txt') }

        def resolved_background = (cleaned_background_ch ?: Channel.empty())
            .collect()
            .ifEmpty {
                def bg = params.scoring_background_input ?: ''
                if (bg && file(bg).exists()) [file(bg)] else [file('NO_BACKGROUND')]
            }
            .map { it && it.size() > 0 ? it[0] : file('NO_BACKGROUND') }

        // ── Run scoring — single pass on full postproc pool ────────────────
        def compute_out = SCORING_COMPUTE(
            resolved_postproc,
            resolved_fade_top,
            resolved_fade_bottom,
            resolved_fade_site_top_ch,
            resolved_fade_site_bot_ch,
            resolved_rer,
            accum_all_ch
        )

        // ── Render report ──────────────────────────────────────────────────
        // Helper: for optional compute outputs, emit a sentinel when absent.
        def _opt = { ch, sentinel ->
            ch.ifEmpty { file(sentinel) }
        }

        // CAAS permulation null (corStat_byrank rds) — resolved once, shared by the
        // scoring report (permulation-null overview) and the FCS report (p.perm).
        // Value channel via collect()/map (NOT .first(), which warns on a value channel).
        //
        // Three ways to get here, in precedence order:
        //   1. caas_perms_ch      — CAAS_PERMULATION ran live in this pass; always valid.
        //   2. caas_pos_detail_file — no live null, but a prior run's raw per-cycle
        //      detail is available: REBUILD rather than import a cached RDS. The null
        //      must hold the same gene-level statistic as the observed gene score
        //      (see scoring_compute.R section 4a); a cached caas_perms.rds carries no
        //      such guarantee, and rebuilding costs minutes because no ASR replay is
        //      involved. fcs_enrich.R detects the mismatch and NAs p.perm otherwise,
        //      so this is the difference between working p.perm and no p.perm.
        //   3. caas_perms_file    — import as-is; only safe if it was produced by the
        //      current scoring formula.
        def caas_perms_resolved
        if (!caas_perms_ch && params.caas_pos_detail_file) {
            def _detail = file(params.caas_pos_detail_file)
            assert _detail.exists() : "SCORING: --caas_pos_detail_file not found: ${params.caas_pos_detail_file}"
            log.info "SCORING: rebuilding CAAS permulation null from ${_detail.name} (no ASR replay)"
            caas_perms_resolved = CAAS_PERMS_REBUILD(Channel.value(_detail), resolved_background)
                .perms
                .collect()
                .map { it[0] }
        } else {
            caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
                .ifEmpty { file(params.caas_perms_file ?: 'NO_FILE') }
                .collect()
                .map { it[0] }
        }

        // Lean per-scheme null files → report permulation section.
        // Fall back to --caas_pos_pval_file / --caas_pos_sample_file so a pass that
        // did not run CAAS_PERMULATION itself (no live CT) can still reuse a prior
        // run's position-level null, mirroring how caas_perms_file is resolved above.
        def caas_pos_pval_resolved = (caas_pos_pval_ch ?: Channel.empty())
            .ifEmpty { file(params.caas_pos_pval_file ?: 'NO_CAAS_POS_PVAL') }
            .collect()
            .map { it[0] }
        def caas_pos_sample_resolved = (caas_pos_sample_ch ?: Channel.empty())
            .ifEmpty { file(params.caas_pos_sample_file ?: 'NO_CAAS_POS_SAMPLE') }
            .collect()
            .map { it[0] }
        def caas_pos_quantiles_resolved = (caas_pos_quantiles_ch ?: Channel.empty())
            .ifEmpty { file(params.caas_pos_quantiles_file ?: 'NO_CAAS_POS_QUANTILES') }
            .collect()
            .map { it[0] }

        // filtered_discovery (resolved_postproc, a value channel) = the observed
        // per-(gene,position,scheme) asr_path_score overlaid on the per-scheme null;
        // when scoring runs it is always present (same condition as position_scores).
        def report_out = SCORING_REPORT(
            compute_out.position_scores,
            compute_out.gene_scores,
            compute_out.gene_correlations,
            _opt(compute_out.stress_summary,         'NO_SCORING_STRESS_SUMMARY'),
            _opt(compute_out.stress_correlations,    'NO_SCORING_STRESS_CORR'),
            _opt(compute_out.stress_rank_agreement,  'NO_SCORING_STRESS_RANK'),
            _opt(compute_out.stress_top_overlap,     'NO_SCORING_STRESS_OVERLAP'),
            _opt(compute_out.stress_variants,        'NO_SCORING_STRESS_VARIANTS'),
            resolved_fade_site_top_ch,
            resolved_fade_site_bot_ch,
            resolved_genomic_info,
            caas_perms_resolved,
            caas_pos_pval_resolved,
            caas_pos_sample_resolved,
            caas_pos_quantiles_resolved,
            resolved_postproc,
            resolved_background
        )

        def final_reports = report_out.report

    emit:
        position_scores  = compute_out.position_scores
        gene_scores      = compute_out.gene_scores
        fcs_stats        = compute_out.fcs_stats
        fcs_stats_rer    = compute_out.fcs_stats_rer
        fcs_stats_fade   = compute_out.fcs_stats_fade
        fcs_stats_accum  = compute_out.fcs_stats_accum
        gene_lists       = compute_out.gene_lists
        position_lists   = compute_out.position_lists
        report           = final_reports
        // The CAAS null actually used here — which may have been REBUILT from
        // --caas_pos_detail_file rather than imported. ENRICHMENT must consume this
        // one, not re-resolve from params, or its FCS p.perm would fall back to the
        // cached (possibly stale) caas_perms.rds while the scoring report used the
        // rebuilt one — two different nulls in the same run.
        caas_perms       = caas_perms_resolved
}
