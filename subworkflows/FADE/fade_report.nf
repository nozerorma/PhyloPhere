#!/usr/bin/env nextflow

/*
 * FADE_REPORT
 * ───────────
 * Generate an HTML summary report from all FADE JSON results for a given
 * direction (top / bottom). Calls FADE_report.Rmd.
 *
 * Inputs
 * ──────
 *   direction  : val — 'top' or 'bottom'
 *   json_files : collected list of *.FADE.json paths
 *
 * Outputs
 * ───────
 *   report      : HTML report
 *   summary_tsv : gene-level summary table (TSV)
 */

process FADE_REPORT {
    tag "fade_report|${direction}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_summary_*.tsv'

    input:
    val  direction
    path json_files

    output:
    path "FADE_report_${direction}.html", emit: report
    path "fade_summary_${direction}.tsv", emit: summary_tsv, optional: true

    script:
    def local_dir   = "${baseDir}/subworkflows/FADE/local"
    def outdir      = "${params.outdir}/selection/fade/${direction}"
    def bf_thr      = params.fade_bf_threshold          ?: 100
    def min_genes   = params.fade_min_genes_for_heatmap ?: 3
    def traitname   = params.traitname ?: 'unknown_trait'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'FADE_report.Rmd',
                params = list(
                    json_dir        = '.',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}'
                ),
                output_file = 'FADE_report_${direction}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'FADE_report.Rmd',
                params = list(
                    json_dir        = '.',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}'
                ),
                output_file = 'FADE_report_${direction}.html'
            )
        "
        """
    }
}
