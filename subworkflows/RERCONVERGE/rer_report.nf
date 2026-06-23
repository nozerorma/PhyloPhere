#!/usr/bin/env nextflow

/*
 * RER_REPORT
 * ──────────
 * Generate an HTML summary report from the RERconverge continuous analysis
 * RDS result file. Calls RERconverge_report.Rmd.
 *
 * Inputs
 * ──────
 *   continuous_output : path — {traitname}.continuous.output RDS file
 *   gmt_file          : path — GMT gene-set annotation file (or 'NO_FILE' to skip enrichment)
 *
 * Outputs
 * ───────
 *   report      : HTML report
 *   summary_tsv : gene-level summary table (TSV, optional)
 */

process RER_REPORT {
    tag "rer_report|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/RERConverge/RER_Results",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/RERConverge/RER_Results",
               mode: 'copy', overwrite: true,
               pattern: 'rerconverge_summary_*.tsv'

    input:
    path continuous_output

    output:
    path "RERconverge_report.html",   emit: report
    path "rerconverge_summary_*.tsv", emit: summary_tsv, optional: true

    script:
    def local_dir        = "${baseDir}/subworkflows/RERCONVERGE/local"
    def outdir           = "${params.outdir}/RERConverge/RER_Results"
    def pval_thr         = params.rer_pval_threshold  ?: 0.05
    def top_n            = params.rer_top_n_labels    ?: 15
    def traitname        = params.traitname           ?: 'unknown_trait'
    def top_pct          = params.rer_top_pct          ?: 0.10
    def top5_pct         = params.rer_top5_pct         ?: 0.05
    def top1_pct         = params.rer_top1_pct         ?: 0.01

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        echo "[RER_REPORT] Input RDS: ${continuous_output}"

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    top_n_labels      = ${top_n},
                    top_pct           = ${top_pct},
                    top5_pct          = ${top5_pct},
                    top1_pct          = ${top1_pct},
                    output_dir        = '${outdir}'
                ),
                output_file = 'RERconverge_report.html'
            )
        "

        if [ -f 'RERconverge_report.html' ]; then
            echo "[RER_REPORT] Report generated: RERconverge_report.html"
        else
            echo "[RER_REPORT] WARNING: Report file was not created."
        fi
        """
    } else {
        """
        cp -R ${local_dir}/* .

        echo "[RER_REPORT] Input RDS: ${continuous_output}"

        Rscript -e "
            rmarkdown::render(
                'RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    top_n_labels      = ${top_n},
                    top_pct           = ${top_pct},
                    top5_pct          = ${top5_pct},
                    top1_pct          = ${top1_pct},
                    output_dir        = '${outdir}'
                ),
                output_file = 'RERconverge_report.html'
            )
        "

        if [ -f 'RERconverge_report.html' ]; then
            echo "[RER_REPORT] Report generated: RERconverge_report.html"
        else
            echo "[RER_REPORT] WARNING: Report file was not created."
        fi
        """
    }
}
