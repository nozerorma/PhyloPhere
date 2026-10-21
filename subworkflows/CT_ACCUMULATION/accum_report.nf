#!/usr/bin/env nextflow

// accum_report.nf — Render the HTML accumulation report.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_REPORT: renders 10.Accumulation_report.Rmd from the accumulation results.
 *
 *  The results CSVs arrive staged flat in the task directory, where the Rmd reads them. The two
 *  script branches differ only in the container entrypoint.
 *
 *  Consumes:  all accumulation_<direction>_results.csv files (collected)
 *  Produces:  10.Accumulation_report.html, accumulation_summary_*.tsv (optional)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Report rendering ─────────────────────────────────────────────────────────

process ACCUMULATION_REPORT {
    tag "accumulation_report|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/accumulation",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/accumulation",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path results_csvs   // all accumulation_<direction>_results.csv files staged flat

    output:
    path "10.Accumulation_report.html",      emit: report
    path "accumulation_summary_*.tsv",    emit: summary_tsv, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def outdir    = "${params.outdir}/accumulation"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr_thr   = params.accumulation_fdr ?: 0.1
    def wrapper   = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh ' : ''
    """
    cp -R ${local_dir}/* .
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

    ${wrapper}Rscript -e "
        rmarkdown::render(
            '10.Accumulation_report.Rmd',
            params = list(
                accum_dir      = '.',
                traitname      = '${traitname}',
                fdr_threshold  = ${fdr_thr},
                output_dir     = '${outdir}'
            ),
            output_file = '10.Accumulation_report.html'
        )
    "
    """
}
