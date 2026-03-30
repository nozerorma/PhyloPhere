#!/usr/bin/env nextflow

/*
 * SCORING_REPORT
 * ──────────────
 * Render an HTML summary report for position-level and gene-level
 * CAAS scores.  Calls SCORING_report.Rmd.
 *
 * Inputs
 * ──────
 *   position_scores  : path — position_scores.tsv
 *   gene_scores      : path — gene_scores.tsv
 *   gene_correlations: path — gene_correlations.tsv
 *
 * Outputs
 * ───────
 *   report : HTML report
 */

process SCORING_REPORT {
    tag "scoring_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path position_scores
    path gene_scores
    path gene_correlations
    path vep_transvar    // optional: TransVar annotation TSV (NO_VEP_TRANSVAR sentinel when absent)
    path vep_primateai   // optional: PrimateAI-3D score TSV  (NO_VEP_PRIMATEAI sentinel when absent)
    path genomic_info    // optional: gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)

    output:
    path "SCORING_report.html", emit: report

    script:
    def local_dir      = "${baseDir}/subworkflows/SCORING/local"
    def outdir         = "${params.outdir}/scoring"
    def traitname      = params.traitname ?: 'unknown_trait'
    def top_pct        = params.scoring_position_top_pct  ?: 0.10
    def top5_pct       = params.scoring_position_top5_pct ?: 0.05
    def top1_pct       = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct      = params.scoring_gene_top_pct      ?: 0.10
    def g_top5_pct     = params.scoring_gene_top5_pct     ?: 0.05
    def g_top1_pct     = params.scoring_gene_top1_pct     ?: 0.01
    // Resolve optional sentinel files: pass 'NULL' (R NULL) when no real file is staged
    def tv_arg  = (vep_transvar.name  =~ /^NO_VEP_TRANSVAR/)  ? 'NULL' : "'${vep_transvar}'"
    def pai_arg = (vep_primateai.name =~ /^NO_VEP_PRIMATEAI/) ? 'NULL' : "'${vep_primateai}'"
    def gi_arg  = (genomic_info.name  =~ /^NO_GENOMIC_INFO/)  ? 'NULL' : "'${genomic_info}'"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'SCORING_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    top5_pct             = ${top5_pct},
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top5_pct        = ${g_top5_pct},
                    gene_top1_pct        = ${g_top1_pct},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    genomic_info_file    = ${gi_arg}
                ),
                output_file = 'SCORING_report.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'SCORING_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    top5_pct             = ${top5_pct},
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top5_pct        = ${g_top5_pct},
                    gene_top1_pct        = ${g_top1_pct},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    genomic_info_file    = ${gi_arg}
                ),
                output_file = 'SCORING_report.html'
            )
        "
        """
    }
}
