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
    path gmt_file           // pass file('NO_FILE') when not provided

    output:
    path "RERconverge_report.html",   emit: report
    path "rerconverge_summary_*.tsv", emit: summary_tsv, optional: true

    script:
    def local_dir        = "${baseDir}/subworkflows/RERCONVERGE/local"
    def outdir           = "${params.outdir}/RERConverge/RER_Results"
    def pval_thr         = params.rer_pval_threshold  ?: 0.05
    def rho_thr          = params.rer_rho_threshold   ?: 0.3
    def top_n            = params.rer_top_n_labels    ?: 15
    def traitname        = params.traitname           ?: 'unknown_trait'
    def enrich_min_genes = params.rer_enrich_min_genes ?: 10
    def enrich_top_n     = params.rer_enrich_top_n     ?: 20
    // gmt_file is staged by Nextflow — pass basename if real, empty string if sentinel
    def gmt_arg          = (gmt_file.name == 'NO_FILE') ? '' : gmt_file.name

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        echo "[RER_REPORT] Input RDS: ${continuous_output}"
        echo "[RER_REPORT] GMT file: ${gmt_arg ?: 'none'}"

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    rho_threshold     = ${rho_thr},
                    top_n_labels      = ${top_n},
                    output_dir        = '${outdir}',
                    gmt_file          = '${gmt_arg}',
                    enrich_min_genes  = ${enrich_min_genes},
                    enrich_top_n      = ${enrich_top_n}
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
        echo "[RER_REPORT] GMT file: ${gmt_arg ?: 'none'}"

        Rscript -e "
            rmarkdown::render(
                'RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    rho_threshold     = ${rho_thr},
                    top_n_labels      = ${top_n},
                    output_dir        = '${outdir}',
                    gmt_file          = '${gmt_arg}',
                    enrich_min_genes  = ${enrich_min_genes},
                    enrich_top_n      = ${enrich_top_n}
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
