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
 *   rer_summary           : path — rerconverge_summary_{trait}.tsv (or NO_FILE)
 *   accum_files           : path — directory or collected CSVs (or NO_FILE)
 *
 * Outputs
 * ───────
 *   position_scores    : TSV with per Gene×Position scores
 *   gene_scores        : TSV with per Gene scores + directional significance flags
 *   gene_correlations  : TSV with pairwise correlations
 *   position_gene_lists: TXT files for position-level ORA (all/top/bottom × 10/5/1%)
 *   gene_gene_lists    : TXT files for gene-level ORA (significance sets + intersections)
 */

process SCORING_COMPUTE {
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}"
    label 'error_retry'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/scoring/gene_lists/position",
               mode: 'copy', overwrite: true,
               pattern: 'pos_gene_lists/*.txt'
    publishDir path: "${params.outdir}/scoring/gene_lists/gene",
               mode: 'copy', overwrite: true,
               pattern: 'gene_gene_lists/*.txt'

    input:
    path postproc_file
    path fade_summary_top
    path fade_summary_bottom
    path rer_summary
    path accum_files

    output:
    path "position_scores.tsv",                              emit: position_scores
    path "gene_scores.tsv",                                  emit: gene_scores
    path "gene_correlations.tsv",                            emit: gene_correlations
    path "position_score_stress_summary.tsv",         optional: true, emit: stress_summary
    path "position_score_stress_correlations.tsv",    optional: true, emit: stress_correlations
    path "position_score_stress_rank_agreement.tsv",  optional: true, emit: stress_rank_agreement
    path "position_score_stress_top_overlap.tsv",     optional: true, emit: stress_top_overlap
    path "position_score_stress_variants.tsv",        optional: true, emit: stress_variants
    path "position_score_stress_latent_loadings.tsv", optional: true, emit: stress_latent_loadings
    path "pos_gene_lists/*.txt",                      optional: true, emit: position_gene_lists
    path "gene_gene_lists/*.txt",                     optional: true, emit: gene_gene_lists

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local"
    def top_pct    = params.scoring_position_top_pct  ?: 0.10
    def top5_pct   = params.scoring_position_top5_pct ?: 0.05
    def top1_pct   = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct  = params.scoring_gene_top_pct      ?: 0.10
    def g_top5_pct = params.scoring_gene_top5_pct     ?: 0.05
    def g_top1_pct = params.scoring_gene_top1_pct     ?: 0.01
    def accum_arg  = (accum_files instanceof List
                        ? (accum_files.size() == 1 && accum_files[0].name.startsWith('NO_') ? 'NO_ACCUM' : '.')
                        : (accum_files.name.startsWith('NO_') ? 'NO_ACCUM' : '.'))

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/scoring_compute.R .

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc    '${postproc_file}' \
            --fade_top    '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer         '${rer_summary}' \
            --accum_dir   '${accum_arg}' \
            --stress      '${params.scoring_stress ?: false}' \
            --stress_top_n ${params.scoring_stress_top_n ?: 25} \
            --top_pct     ${top_pct} \
            --top5_pct    ${top5_pct} \
            --top1_pct    ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top5_pct ${g_top5_pct} \
            --gene_top1_pct ${g_top1_pct}
        """
    } else {
        """
        cp ${local_dir}/scoring_compute.R .

        Rscript scoring_compute.R \
            --postproc    '${postproc_file}' \
            --fade_top    '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer         '${rer_summary}' \
            --accum_dir   '${accum_arg}' \
            --stress      '${params.scoring_stress ?: false}' \
            --stress_top_n ${params.scoring_stress_top_n ?: 25} \
            --top_pct     ${top_pct} \
            --top5_pct    ${top5_pct} \
            --top1_pct    ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top5_pct ${g_top5_pct} \
            --gene_top1_pct ${g_top1_pct}
        """
    }
}
