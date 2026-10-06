#!/usr/bin/env nextflow
// rer_report.nf — HTML report and gene-level summary of a RERconverge correlation table.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_REPORT: renders 5.RERconverge_report.Rmd over the RDS written by RER_CONT or
 *  RER_BIN (a table with Rho, N, P and p.adj per gene, and p.perm with permutations).
 *  A failed task is ignored.
 *
 *  Consumes:  correlation RDS (<trait>.continuous.output or <trait>.binary.output)
 *  Produces:  5.RERconverge_report.html, rerconverge_summary_<trait>.tsv (gene-level table,
 *             emitted when the report writes it)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Report ─────────────────────────────────────────────────────────────────────

process RER_REPORT {
    tag "rer_report|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/rerconverge/rer_results",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/rerconverge/rer_results",
               mode: 'copy', overwrite: true,
               pattern: 'rerconverge_summary_*.tsv'

    input:
    path continuous_output

    output:
    path "5.RERconverge_report.html",   emit: report
    path "rerconverge_summary_*.tsv", emit: summary_tsv, optional: true

    script:
    def local_dir        = "${baseDir}/subworkflows/RERCONVERGE/local"
    def outdir           = "${params.outdir}/rerconverge/rer_results"
    def pval_thr         = params.rer_pval_threshold  ?: 0.05
    def top_n            = params.rer_top_n_labels    ?: 15
    def traitname        = params.traitname           ?: 'unknown_trait'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        echo "[RER_REPORT] Input RDS: ${continuous_output}"

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '5.RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    top_n_labels      = ${top_n},
                    output_dir        = '${outdir}'
                ),
                output_file = '5.RERconverge_report.html'
            )
        "

        if [ -f '5.RERconverge_report.html' ]; then
            echo "[RER_REPORT] Report generated: 5.RERconverge_report.html"
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
                '5.RERconverge_report.Rmd',
                params = list(
                    continuous_rds    = '${continuous_output}',
                    traitname         = '${traitname}',
                    pval_threshold    = ${pval_thr},
                    top_n_labels      = ${top_n},
                    output_dir        = '${outdir}'
                ),
                output_file = '5.RERconverge_report.html'
            )
        "

        if [ -f '5.RERconverge_report.html' ]; then
            echo "[RER_REPORT] Report generated: 5.RERconverge_report.html"
        else
            echo "[RER_REPORT] WARNING: Report file was not created."
        fi
        """
    }
}
