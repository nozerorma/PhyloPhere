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

    output:
    path "SCORING_report.html", emit: report

    script:
    def local_dir = "${baseDir}/subworkflows/SCORING/local"
    def outdir    = "${params.outdir}/scoring"
    def traitname = params.traitname ?: 'unknown_trait'
    def top_pct   = params.scoring_position_top_pct  ?: 0.10
    def top1_pct  = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct = params.scoring_gene_top_pct      ?: 0.10
    def g_top1_pct = params.scoring_gene_top1_pct    ?: 0.01

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
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top1_pct        = ${g_top1_pct}
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
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top1_pct        = ${g_top1_pct}
                ),
                output_file = 'SCORING_report.html'
            )
        "
        """
    }
}
