#!/usr/bin/env nextflow
// fade_report.nf — HTML report of the FADE results of one direction.
// PhyloPhere | subworkflows/FADE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  FADE_REPORT: renders 6.FADE_report.Rmd over the merged FADE tables of one direction
 *  (top or bottom), as written by FADE_TABLES.
 *
 *  Consumes:  direction, fade_gene_counts_<direction>.tsv, fade_sites_<direction>.csv,
 *             optional foreground species list (NO_FG_LIST sentinel when absent)
 *  Produces:  6.FADE_report_<direction>.html (optional: a failed render does not stop the run)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Report ─────────────────────────────────────────────────────────────────────

process FADE_REPORT {
    tag "fade_report|${direction}"
    label 'process_reporting'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    val  direction
    path gene_counts_file
    path sites_file
    path fg_list_file  // optional: foreground species list (NO_FG_LIST sentinel when absent)

    output:
    path "6.FADE_report_${direction}.html", emit: report, optional: true

    script:
    def local_dir   = "${baseDir}/subworkflows/FADE/local"
    def outdir      = "${params.outdir}/selection/fade/${direction}"
    def bf_thr      = params.fade_bf_threshold          ?: 100
    def min_genes   = params.fade_min_genes_for_heatmap ?: 3
    def traitname   = params.traitname ?: 'unknown_trait'
    def fg_arg      = (fg_list_file.name =~ /^NO_FG_LIST/) ? 'NULL' : "'${fg_list_file}'"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        echo "[FADE_REPORT] direction=${direction} | genes in the table: \$(( \$(wc -l < '${gene_counts_file}') - 1 ))"

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '6.FADE_report.Rmd',
                params = list(
                    gene_counts_file = '${gene_counts_file}',
                    sites_file       = '${sites_file}',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}',
                    fg_list_file    = ${fg_arg}
                ),
                output_file = '6.FADE_report_${direction}.html'
            )
        " || echo "[FADE_REPORT] WARNING: rendering stopped early."

        if [ -f '6.FADE_report_${direction}.html' ]; then
            echo "[FADE_REPORT] Report generated: 6.FADE_report_${direction}.html"
        else
            echo "[FADE_REPORT] WARNING: Report file was not created for direction '${direction}'."
        fi
        """
    } else {
        """
        cp -R ${local_dir}/* .

        echo "[FADE_REPORT] direction=${direction} | genes in the table: \$(( \$(wc -l < '${gene_counts_file}') - 1 ))"

        Rscript -e "
            rmarkdown::render(
                '6.FADE_report.Rmd',
                params = list(
                    gene_counts_file = '${gene_counts_file}',
                    sites_file       = '${sites_file}',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}',
                    fg_list_file    = ${fg_arg}
                ),
                output_file = '6.FADE_report_${direction}.html'
            )
        " || echo "[FADE_REPORT] WARNING: rendering stopped early."

        if [ -f '6.FADE_report_${direction}.html' ]; then
            echo "[FADE_REPORT] Report generated: 6.FADE_report_${direction}.html"
        else
            echo "[FADE_REPORT] WARNING: Report file was not created for direction '${direction}'."
        fi
        """
    }
}
