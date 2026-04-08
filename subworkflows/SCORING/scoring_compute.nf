#!/usr/bin/env nextflow

/*
 * SCORING_COMPUTE
 * ───────────────
 * Compute composite CAAS scores at position-level and gene-level.
 * Integrates outputs from CT_POSTPROC, PGLS, FADE, RERConverge,
 * CT_ACCUMULATION, and MoleRate.
 *
 * Inputs
 * ──────
 *   postproc_file         : path — filtered_discovery.tsv (mandatory)
 *   pgls_file             : path — site_pgls.tsv (or NO_FILE sentinel)
 *   fade_summary_top      : path — fade_summary_top.tsv (or NO_FILE)
 *   fade_summary_bottom   : path — fade_summary_bottom.tsv (or NO_FILE)
 *   rer_summary           : path — rerconverge_summary_{trait}.tsv (or NO_FILE)
 *   accum_files           : path — directory or collected CSVs (or NO_FILE)
 *   molerate_summary_top  : path — molerate_summary_top.tsv (or NO_FILE)
 *   molerate_summary_bottom: path — molerate_summary_bottom.tsv (or NO_FILE)
 *
 * Outputs
 * ───────
 *   position_scores    : TSV with per Gene×Position scores
 *   gene_scores        : TSV with per Gene scores
 *   gene_correlations  : TSV with pairwise correlations
 *   position_gene_lists: TXT files for position-level ORA
 *   gene_gene_lists    : TXT files for gene-level ORA
 */

process SCORING_COMPUTE {
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}|${direction}"
    label 'error_retry'

    publishDir path: { "${params.outdir}/scoring/${direction}" },
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: { "${params.outdir}/scoring/${direction}/gene_lists/position" },
               mode: 'copy', overwrite: true,
               pattern: 'pos_gene_lists/*.txt'
    publishDir path: { "${params.outdir}/scoring/${direction}/gene_lists/gene" },
               mode: 'copy', overwrite: true,
               pattern: 'gene_gene_lists/*.txt'

    input:
    val  direction
    path postproc_file
    path pgls_file
    path pgls_excess_file
    path fade_summary_top
    path fade_summary_bottom
    path rer_summary
    path accum_files
    path molerate_summary_top
    path molerate_summary_bottom

    output:
    val  direction,                                             emit: direction
    tuple val(direction), path("position_scores.tsv"),         emit: position_scores
    tuple val(direction), path("gene_scores.tsv"),             emit: gene_scores
    tuple val(direction), path("gene_correlations.tsv"),       emit: gene_correlations
    tuple val(direction), path("pgls_excess_for_report.tsv"),  emit: pgls_excess_report,     optional: true
    tuple val(direction), path("position_score_stress_summary.tsv"),         emit: stress_summary,       optional: true
    tuple val(direction), path("position_score_stress_correlations.tsv"),    emit: stress_correlations,  optional: true
    tuple val(direction), path("position_score_stress_rank_agreement.tsv"),  emit: stress_rank_agreement, optional: true
    tuple val(direction), path("position_score_stress_top_overlap.tsv"),     emit: stress_top_overlap,   optional: true
    tuple val(direction), path("position_score_stress_variants.tsv"),        emit: stress_variants,      optional: true
    tuple val(direction), path("position_score_stress_latent_loadings.tsv"), emit: stress_latent_loadings, optional: true
    path "pos_gene_lists/*.txt",                                emit: position_gene_lists, optional: true
    path "gene_gene_lists/*.txt",                               emit: gene_gene_lists,     optional: true

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local"
    def top_pct    = params.scoring_position_top_pct  ?: 0.10
    def top5_pct   = params.scoring_position_top5_pct ?: 0.05
    def top1_pct   = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct  = params.scoring_gene_top_pct      ?: 0.10
    def g_top5_pct = params.scoring_gene_top5_pct     ?: 0.05
    def g_top1_pct = params.scoring_gene_top1_pct     ?: 0.01
    // Determine accum_dir: if accum_files is a real directory, pass it;
    // otherwise check if multiple CSVs were staged (Nextflow stages them flat).
    // accum_files may be a single sentinel file or a list of CSVs
    def accum_arg  = (accum_files instanceof List
                        ? (accum_files.size() == 1 && accum_files[0].name.startsWith('NO_') ? 'NO_ACCUM' : '.')
                        : (accum_files.name.startsWith('NO_') ? 'NO_ACCUM' : '.'))

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/scoring_compute.R .

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc   '${postproc_file}' \
            --pgls       '${pgls_file}' \
            --pgls_excess '${pgls_excess_file}' \
            --fade_top   '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer        '${rer_summary}' \
            --accum_dir  '${accum_arg}' \
            --molerate_top    '${molerate_summary_top}' \
            --molerate_bottom '${molerate_summary_bottom}' \
            --stress     '${params.scoring_stress ?: false}' \
            --stress_top_n ${params.scoring_stress_top_n ?: 25} \
            --top_pct    ${top_pct} \
            --top5_pct   ${top5_pct} \
            --top1_pct   ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top5_pct ${g_top5_pct} \
            --gene_top1_pct ${g_top1_pct} \
            --direction     '${direction}'

        for f in pos_gene_lists/*.txt gene_gene_lists/*.txt; do
            [ -f "\$f" ] && mv "\$f" "\${f%.txt}_${direction}.txt"
        done
        """
    } else {
        """
        cp ${local_dir}/scoring_compute.R .

        Rscript scoring_compute.R \
            --postproc   '${postproc_file}' \
            --pgls       '${pgls_file}' \
            --pgls_excess '${pgls_excess_file}' \
            --fade_top   '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer        '${rer_summary}' \
            --accum_dir  '${accum_arg}' \
            --molerate_top    '${molerate_summary_top}' \
            --molerate_bottom '${molerate_summary_bottom}' \
            --stress     '${params.scoring_stress ?: false}' \
            --stress_top_n ${params.scoring_stress_top_n ?: 25} \
            --top_pct    ${top_pct} \
            --top5_pct   ${top5_pct} \
            --top1_pct   ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top5_pct ${g_top5_pct} \
            --gene_top1_pct ${g_top1_pct} \
            --direction     '${direction}'

        for f in pos_gene_lists/*.txt gene_gene_lists/*.txt; do
            [ -f "\$f" ] && mv "\$f" "\${f%.txt}_${direction}.txt"
        done
        """
    }
}
