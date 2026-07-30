#!/usr/bin/env nextflow

/*
 * SCORING_COMPUTE
 * ───────────────
 * Compute CAAS scores at position-level and gene-level.
 * Integrates outputs from CT_POSTPROC, FADE, RERConverge and CT_ACCUMULATION.
 * Runs once on the full postproc pool; directional characterisation is handled
 * post-scoring via the change_side column.
 *
 * Inputs
 * ──────
 *   postproc_file         : path — filtered_discovery.tsv (mandatory)
 *   fade_summary_top      : path — fade_summary_top.tsv (or NO_FILE)
 *   fade_summary_bottom   : path — fade_summary_bottom.tsv (or NO_FILE)
 *   fade_site_top         : path — fade_site_bf_top.tsv (or NO_FADE_SITE_TOP sentinel)
 *   fade_site_bot         : path — fade_site_bf_bottom.tsv (or NO_FADE_SITE_BOT sentinel)
 *   rer_summary           : path — rerconverge_summary_{trait}.tsv (or NO_FILE)
 *   accum_files           : path — directory or collected CSVs (or NO_FILE)
 *
 * Outputs
 * ───────
 *   position_scores    : TSV with per Gene×Position scores
 *   gene_scores        : TSV with per Gene scores + directional significance flags
 *   gene_correlations  : TSV with pairwise correlations
 *   gene_lists                : directory of 12 slice TSVs (Top/Bottom/Global × 25/10/5/1%) for AMI/reports
 *   position_lists            : directory of 12 slice TSVs (Top/Bottom/Global × 25/10/5/1%), position-level,
 *                               for posenrich/reports (Gene, Position, CAAS_score)
 *   gene_threshold_enrichment : TSV — gene-level enrichment curve across CAAS thresholds × tools
 *   pos_threshold_enrichment  : TSV — position-level FADE enrichment curve across CAAS thresholds
 */

process SCORING_COMPUTE {
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}"
    label 'error_retry'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: 'gene_lists'
    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: 'position_lists'

    input:
    path postproc_file
    path fade_summary_top
    path fade_summary_bottom
    path fade_site_top
    path fade_site_bot
    path rer_summary
    path accum_files

    output:
    path "position_scores.tsv",                              emit: position_scores
    path "gene_scores.tsv",                                  emit: gene_scores
    path "fcs_stats.tsv",                                    emit: fcs_stats
    path "fcs_stats_rer.tsv",                          optional: true, emit: fcs_stats_rer
    path "fcs_stats_fade.tsv",                         optional: true, emit: fcs_stats_fade
    path "fcs_stats_accum.tsv",                        optional: true, emit: fcs_stats_accum
    path "gene_correlations.tsv",                            emit: gene_correlations
    path "position_score_stress_summary.tsv",         optional: true, emit: stress_summary
    path "position_score_stress_correlations.tsv",    optional: true, emit: stress_correlations
    path "position_score_stress_rank_agreement.tsv",  optional: true, emit: stress_rank_agreement
    path "position_score_stress_top_overlap.tsv",     optional: true, emit: stress_top_overlap
    path "position_score_stress_variants.tsv",        optional: true, emit: stress_variants
    path "position_score_stress_latent_loadings.tsv", optional: true, emit: stress_latent_loadings
    path "gene_lists",                                 optional: true, emit: gene_lists
    path "position_lists",                             optional: true, emit: position_lists
    path "gene_threshold_enrichment.tsv",             optional: true, emit: gene_threshold_enrichment
    path "pos_threshold_enrichment.tsv",              optional: true, emit: pos_threshold_enrichment

    script:
    def local_dir       = "${baseDir}/subworkflows/SCORING/local/src"
    def top_pct         = params.scoring_position_top_pct    ?: 0.10
    def g_top_pct       = params.scoring_gene_top_pct        ?: 0.10
    def accum_arg         = (accum_files instanceof List
                              ? (accum_files.size() == 1 && accum_files[0].name.startsWith('NO_') ? 'NO_ACCUM' : '.')
                              : (accum_files.name.startsWith('NO_') ? 'NO_ACCUM' : '.'))
    def fs_top_arg        = fade_site_top.name =~ /^NO_FADE_SITE_TOP/ ? 'NO_FADE_SITE_TOP' : "${fade_site_top}"
    def fs_bot_arg        = fade_site_bot.name =~ /^NO_FADE_SITE_BOT/ ? 'NO_FADE_SITE_BOT' : "${fade_site_bot}"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/scoring_compute.R .

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc       '${postproc_file}' \
            --fade_top       '${fade_summary_top}' \
            --fade_bottom    '${fade_summary_bottom}' \
            --fade_site_top  '${fs_top_arg}' \
            --fade_site_bot  '${fs_bot_arg}' \
            --rer            '${rer_summary}' \
            --accum_dir      '${accum_arg}' \
            --stress               '${params.scoring_stress ?: false}' \
            --stress_top_n         ${params.scoring_stress_top_n ?: 25} \
            --top_pct              ${top_pct} \
            --top25_pct            0.25 \
            --top5_pct             0.05 \
            --top1_pct             0.01 \
            --gene_top_pct         ${g_top_pct} \
            --gene_top25_pct       0.25 \
            --gene_top5_pct        0.05 \
            --gene_top1_pct        0.01
        """
    } else {
        """
        cp ${local_dir}/scoring_compute.R .

        Rscript scoring_compute.R \
            --postproc       '${postproc_file}' \
            --fade_top       '${fade_summary_top}' \
            --fade_bottom    '${fade_summary_bottom}' \
            --fade_site_top  '${fs_top_arg}' \
            --fade_site_bot  '${fs_bot_arg}' \
            --rer            '${rer_summary}' \
            --accum_dir      '${accum_arg}' \
            --stress               '${params.scoring_stress ?: false}' \
            --stress_top_n         ${params.scoring_stress_top_n ?: 25} \
            --top_pct              ${top_pct} \
            --top25_pct            0.25 \
            --top5_pct             0.05 \
            --top1_pct             0.01 \
            --gene_top_pct         ${g_top_pct} \
            --gene_top25_pct       0.25 \
            --gene_top5_pct        0.05 \
            --gene_top1_pct        0.01
        """
    }
}
