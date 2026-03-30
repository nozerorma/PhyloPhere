#!/usr/bin/env nextflow

/*
 * MOLERATE_REPORT
 * ───────────────
 * Generate an HTML summary report from all MoleRate JSON results for a given
 * direction (top / bottom). Calls MoleRate_report.Rmd.
 *
 * Inputs
 * ──────
 *   direction  : val — 'top' or 'bottom'
 *   json_files : collected list of *.molerate.json paths
 *
 * Outputs
 * ───────
 *   report      : HTML report
 *   summary_tsv : gene-level summary table (TSV)
 */

process MOLERATE_REPORT {
    tag "molerate_report|${direction}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/selection/molerate/${direction}",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/selection/molerate/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'molerate_summary_*.tsv'

    input:
    val  direction
    path json_files

    output:
    path "MoleRate_report_${direction}.html", emit: report
    path "molerate_summary_${direction}.tsv", emit: summary_tsv, optional: true

    script:
    def local_dir       = "${baseDir}/subworkflows/MOLERATE/local"
    def outdir          = "${params.outdir}/selection/molerate/${direction}"
    def pval_thr        = params.molerate_pval_threshold  ?: 0.05
    def log2fc_thr      = params.molerate_log2fc_threshold ?: 0.5
    def top_n           = params.molerate_top_n_labels    ?: 15
    def traitname       = params.traitname ?: 'unknown_trait'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'MoleRate_report.Rmd',
                params = list(
                    json_dir      = '.',
                    direction     = '${direction}',
                    traitname     = '${traitname}',
                    pval_threshold  = ${pval_thr},
                    log2fc_threshold = ${log2fc_thr},
                    top_n_labels  = ${top_n},
                    output_dir    = '${outdir}'
                ),
                output_file = 'MoleRate_report_${direction}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'MoleRate_report.Rmd',
                params = list(
                    json_dir      = '.',
                    direction     = '${direction}',
                    traitname     = '${traitname}',
                    pval_threshold  = ${pval_thr},
                    log2fc_threshold = ${log2fc_thr},
                    top_n_labels  = ${top_n},
                    output_dir    = '${outdir}'
                ),
                output_file = 'MoleRate_report_${direction}.html'
            )
        "
        """
    }
}
