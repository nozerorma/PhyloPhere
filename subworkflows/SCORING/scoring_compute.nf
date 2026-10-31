#!/usr/bin/env nextflow
include { scoringQuote } from './scoring_contract'
// scoring_compute.nf — Position-level and gene-level CAAS scores, joined with FADE and RER evidence.
// PhyloPhere | subworkflows/SCORING/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_COMPUTE: scores the observed positions and genes of the filtered discovery
 *  and integrates the other tools' evidence per gene.
 *
 *  observed_core_scores.py first computes the position and gene CAAS scores with the
 *  code that also scores the permulation null (src/core/scores.py, copied from
 *  CT_DISAMBIGUATION); scoring_compute.R then adds the position and gene permulation
 *  p-values, joins FADE and RERConverge, and writes the tables, the
 *  ranked slices and the enrichment curves. It runs once on the full postproc pool;
 *  direction is carried by the side column.
 *
 *  Consumes:  postproc_file (filtered_discovery.tsv, mandatory); optional FADE gene
 *             summaries and per-site tables, RERConverge summary,
 *             hypothesis pairs, and the counters of the CAAS permulation null
 *             (null_counters.tsv, null_fit.tsv and null_universe_base.tsv of
 *             SCORING_NULL_POSITIONS) with its cycle roster (N). Each optional input
 *             is a staged NO_* sentinel file when absent.
 *  Produces:  position_scores.tsv, gene_scores.tsv, gene_correlations.tsv,
 *             fcs_stats.tsv (plus fcs_stats_{rer,fade}.tsv when that evidence
 *             exists), null_universe_summary.tsv (N, size of the positional universe,
 *             observed positions the null never scored, saturation curve of the
 *             null; written when a null is given). The nominal activity lists
 *             (gene_lists/, position_lists/) are produced by SCORING_EVIDENCE.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Scoring ────────────────────────────────────────────────────────────────────

process SCORING_COMPUTE {
    cache 'deep'
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}"
    label 'error_retry'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv',
               saveAs: { name -> name in ['gene_scores.tsv','position_scores.tsv','fcs_stats.tsv'] ? null : name }

    input:
    val source_hash
    path postproc_file
    path fade_summary_top
    path fade_summary_bottom
    path fade_site_top
    path fade_site_bot
    path rer_summary
    path hypotheses_pairs   // optional: contrast_hypotheses_pairs.tsv (FOP); NO_HYP_PAIRS sentinel otherwise
    path null_counters      // optional: null_counters.tsv, the counters of the position-level null (SCORING_NULL_POSITIONS); NO_* sentinel otherwise
    path null_fit           // optional: null_fit.tsv, legacy staging metadata (unused); NO_* sentinel otherwise
    path null_base          // optional: null_universe_base.tsv, what the null fixes of the universe summary; NO_* sentinel otherwise
    path cycle_roster       // optional: cycle_roster.txt, N of the null (cycles with no detection included); NO_* sentinel otherwise

    output:
    path "position_scores.tsv",                              emit: position_scores
    path "gene_scores.tsv",                                  emit: gene_scores
    path "fcs_stats.tsv",                                    emit: fcs_stats
    path "fcs_stats_rer.tsv",                          optional: true, emit: fcs_stats_rer
    path "fcs_stats_fade.tsv",                         optional: true, emit: fcs_stats_fade
    path "gene_correlations.tsv",                            emit: gene_correlations
    path "null_universe_summary.tsv",                 optional: true, emit: null_universe_summary

    script:
    def local_dir       = "${baseDir}/subworkflows/SCORING/local/src"
    def core_dir        = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core"
    def fs_top_arg        = fade_site_top.name =~ /^NO_FADE_SITE_TOP/ ? 'NO_FADE_SITE_TOP' : "${fade_site_top}"
    def fs_bot_arg        = fade_site_bot.name =~ /^NO_FADE_SITE_BOT/ ? 'NO_FADE_SITE_BOT' : "${fade_site_bot}"
    def hp_arg            = hypotheses_pairs.name =~ /^NO_/ ? 'NO_HYP_PAIRS' : "${hypotheses_pairs}"
    def counters_arg      = null_counters.name =~ /^NO_/ ? 'NO_FILE' : "${null_counters}"
    def fit_arg           = null_fit.name =~ /^NO_/ ? 'NO_FILE' : "${null_fit}"
    def base_arg          = null_base.name =~ /^NO_/ ? 'NO_FILE' : "${null_base}"
    def roster_arg        = cycle_roster.name =~ /^NO_/ ? 'NO_FILE' : "${cycle_roster}"
    def gene_stat         = params.gene_stat ?: 'max'

    if (params.use_singularity || params.use_apptainer) {
        """
        export VROOM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
        echo '${source_hash}' > scoring_contract.txt
        cp ${local_dir}/scoring_compute.R ${local_dir}/scoring_io.R ${local_dir}/scoring_input.py ${local_dir}/aa_grouping.R ${local_dir}/null_position_p.R ${local_dir}/observed_core_scores.py .
        mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py

        /usr/local/bin/_entrypoint.sh python3 observed_core_scores.py \
            --input ${scoringQuote(postproc_file)} \
            --positions-out core_positions.tsv \
            --genes-out core_genes.tsv \
            --gene-stat ${gene_stat}

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc       ${scoringQuote(postproc_file)} \
            --core_positions core_positions.tsv \
            --core_genes     core_genes.tsv \
            --fade_top       ${scoringQuote(fade_summary_top)} \
            --fade_bottom    ${scoringQuote(fade_summary_bottom)} \
            --fade_site_top  ${scoringQuote(fs_top_arg)} \
            --fade_site_bot  ${scoringQuote(fs_bot_arg)} \
            --rer            ${scoringQuote(rer_summary)} \
            --hypotheses_pairs ${scoringQuote(hp_arg)} \
            --null_counters   ${scoringQuote(counters_arg)} \
            --null_fit        ${scoringQuote(fit_arg)} \
            --null_base       ${scoringQuote(base_arg)} \
            --cycle_roster    ${scoringQuote(roster_arg)} \
            --nominal_threshold ${params.evidence_p_threshold ?: 0.05} \
            --fdr_threshold ${params.evidence_fdr_threshold ?: 0.05} \
            --fade_bf_threshold ${params.fade_bf_threshold ?: 100}
        """
    } else {
        """
        export VROOM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
        echo '${source_hash}' > scoring_contract.txt
        cp ${local_dir}/scoring_compute.R ${local_dir}/scoring_io.R ${local_dir}/scoring_input.py ${local_dir}/aa_grouping.R ${local_dir}/null_position_p.R ${local_dir}/observed_core_scores.py .
        mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py

        python3 observed_core_scores.py \
            --input ${scoringQuote(postproc_file)} \
            --positions-out core_positions.tsv \
            --genes-out core_genes.tsv \
            --gene-stat ${gene_stat}

        Rscript scoring_compute.R \
            --postproc       ${scoringQuote(postproc_file)} \
            --core_positions core_positions.tsv \
            --core_genes     core_genes.tsv \
            --fade_top       ${scoringQuote(fade_summary_top)} \
            --fade_bottom    ${scoringQuote(fade_summary_bottom)} \
            --fade_site_top  ${scoringQuote(fs_top_arg)} \
            --fade_site_bot  ${scoringQuote(fs_bot_arg)} \
            --rer            ${scoringQuote(rer_summary)} \
            --hypotheses_pairs ${scoringQuote(hp_arg)} \
            --null_counters   ${scoringQuote(counters_arg)} \
            --null_fit        ${scoringQuote(fit_arg)} \
            --null_base       ${scoringQuote(base_arg)} \
            --cycle_roster    ${scoringQuote(roster_arg)} \
            --nominal_threshold ${params.evidence_p_threshold ?: 0.05} \
            --fdr_threshold ${params.evidence_fdr_threshold ?: 0.05} \
            --fade_bf_threshold ${params.fade_bf_threshold ?: 100}
        """
    }
}
