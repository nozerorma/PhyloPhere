#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: CAAS Scoring Workflow
 *
 * Computes composite position-level and gene-level CAAS scores by
 * integrating outputs from CT_POSTPROC, FADE, RERConverge,
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
        fade_summary_top_ch      // Channel<path> or null — fade_summary_top.tsv
        fade_summary_bottom_ch   // Channel<path> or null — fade_summary_bottom.tsv
        rer_summary_ch           // Channel<path> or null — rerconverge_summary_{trait}.tsv
        accum_ch                 // Channel<path> or null — collected accumulation CSVs
        molerate_summary_top_ch  // Channel<path> or null — molerate_summary_top.tsv
        molerate_summary_bot_ch  // Channel<path> or null — molerate_summary_bottom.tsv
        background_ch            // Channel<path> or null — cleaned_background_main.txt
        vep_transvar_ch          // Channel<path> or null — transvar_annotations.tsv (optional)
        vep_primateai_ch         // Channel<path> or null — primateai_scores.tsv     (optional)
        vep_aa2prot_ch           // Channel<path> or null — aa2prot_global.csv (alignment→protein pos map)
        genomic_info_ch          // Channel<path> or null — gene genomic coords TSV  (optional)

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

        // Accumulation: collect all CSVs (top + bottom) into a single value channel.
        // IMPORTANT: do NOT use .combine(accum_ch) in the multiMap chain below —
        // when a channel emits List<Path> (from a glob output), .combine() spreads the
        // list as individual tuple elements, causing an arity mismatch in multiMap.
        // Instead we flatten + collect here so the process input receives all files at once
        // and the R script selects the right direction via --direction / --accum_dir '.'.
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

        def resolved_molerate_top = (molerate_summary_top_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_molerate_summary_top ?: 'NO_MOLERATE_TOP') }

        def resolved_molerate_bot = (molerate_summary_bot_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_molerate_summary_bottom ?: 'NO_MOLERATE_BOTTOM') }

        // VEP optional inputs — sentinels prevent staging collisions
        def resolved_vep_transvar = (vep_transvar_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_transvar ?: 'NO_VEP_TRANSVAR') }

        def resolved_vep_primateai = (vep_primateai_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_primateai ?: 'NO_VEP_PRIMATEAI') }

        def resolved_vep_aa2prot = (vep_aa2prot_ch ?: Channel.empty())
            .ifEmpty { file(params.scoring_vep_aa2prot ?: 'NO_VEP_AA2PROT') }

        // Genomic info (gene coordinates) — reuses params.gene_ensembl_file if available
        def resolved_genomic_info = (genomic_info_ch ?: Channel.empty())
            .ifEmpty {
                def gi = params.gene_ensembl_file ?: ''
                if (gi && file(gi).exists()) file(gi) else file('NO_GENOMIC_INFO')
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

        // ── Run scoring computation for each phenotype direction ──────────
        // Each resolved_* is a single-item queue channel (one file per run).
        // combine() creates the Cartesian product: 2 directions × 1 file = 2 items.
        // multiMap then splits the combined tuples back into named channels,
        // giving SCORING_COMPUTE two properly-paired calls without consuming any
        // channel more than once and without using the deprecated .first() operator.
        def directions_ch = Channel.of("top", "bottom")

        def cmp_in = directions_ch
            .combine(resolved_postproc)
            .combine(resolved_fade_top)
            .combine(resolved_fade_bottom)
            .combine(resolved_rer)
            .combine(resolved_molerate_top)
            .combine(resolved_molerate_bot)
            .multiMap { dir, postproc, fade_t, fade_b, rer, mol_t, mol_b ->
                direction:   dir
                postproc:    postproc
                fade_top:    fade_t
                fade_bot:    fade_b
                rer:         rer
                mol_top:     mol_t
                mol_bot:     mol_b
            }

        def compute_out = SCORING_COMPUTE(
            cmp_in.direction,
            cmp_in.postproc,
            cmp_in.fade_top,
            cmp_in.fade_bot,
            cmp_in.rer,
            accum_all_ch,
            cmp_in.mol_top,
            cmp_in.mol_bot
        )

        // ── Render one report per direction ────────────────────────────────
        // SCORING_COMPUTE emits (direction, file) tuples for all outputs.
        // Join on the direction key so each report receives correctly-paired inputs.
        // For optional outputs, join with remainder:true and map null → sentinel file.

        def report_core = compute_out.position_scores   // (dir, pos)
            .join(compute_out.gene_scores)              // (dir, pos, gene)
            .join(compute_out.gene_correlations)        // (dir, pos, gene, corr)

        // For each direction-keyed optional channel, emit a sentinel when absent.
        // report_core.map extracts just the direction keys as the left side of join.
        def _opt = { ch, sentinel ->
            report_core.map { d, _p, _g, _c -> d }
                .join(ch, remainder: true)
                .map { dir, f -> [dir, f ?: file(sentinel)] }
        }

        // Build a fully direction-keyed channel with all report inputs
        def report_all = report_core
            .join(_opt(compute_out.stress_summary,         'NO_SCORING_STRESS_SUMMARY'))
            .join(_opt(compute_out.stress_correlations,    'NO_SCORING_STRESS_CORR'))
            .join(_opt(compute_out.stress_rank_agreement,  'NO_SCORING_STRESS_RANK'))
            .join(_opt(compute_out.stress_top_overlap,     'NO_SCORING_STRESS_OVERLAP'))
            .join(_opt(compute_out.stress_variants,        'NO_SCORING_STRESS_VARIANTS'))
            .join(_opt(compute_out.stress_latent_loadings, 'NO_SCORING_STRESS_LOADINGS'))
        // report_all: (dir, pos, gene, corr, ss, sc, sr, so, sv, sl) — 10 elements

        // Combine with shared VEP/genomic inputs (single-item channels reused per direction)
        def rpt_in = report_all
            .combine(resolved_vep_transvar)
            .combine(resolved_vep_primateai)
            .combine(resolved_vep_aa2prot)
            .combine(resolved_genomic_info)
            .multiMap { dir, pos, gene, corr, ss, sc, sr, so, sv, sl, tv, pai, a2p, gi ->
                direction:   dir
                pos_scores:  pos
                gene_scores: gene
                gene_corr:   corr
                stress_sum:  ss
                stress_corr: sc
                stress_rank: sr
                stress_over: so
                stress_var:  sv
                stress_load: sl
                transvar:    tv
                primateai:   pai
                aa2prot:     a2p
                genomic:     gi
            }

        def report_out = SCORING_REPORT(
            rpt_in.direction,
            rpt_in.pos_scores,
            rpt_in.gene_scores,
            rpt_in.gene_corr,
            rpt_in.stress_sum,
            rpt_in.stress_corr,
            rpt_in.stress_rank,
            rpt_in.stress_over,
            rpt_in.stress_var,
            rpt_in.stress_load,
            rpt_in.transvar,
            rpt_in.primateai,
            rpt_in.aa2prot,
            rpt_in.genomic
        )

        // ── ORA on scoring gene lists (optional) ───────────────────────────
        if (params.scoring_ora) {
            def bg_after_report = resolved_background
                .combine(report_out.report.collect())
                .map { it[0] }  // bg is always first; avoid destructuring (collect() spreads list elements)

            ORA_SCORING_POSITION(
                bg_after_report,
                compute_out.position_gene_lists.collect()
            )

            ORA_SCORING_GENE(
                bg_after_report,
                compute_out.gene_gene_lists.collect()
            )
        }

    emit:
        position_scores = compute_out.position_scores.map { _dir, f -> f }
        gene_scores     = compute_out.gene_scores.map     { _dir, f -> f }
        report          = report_out.report
}
