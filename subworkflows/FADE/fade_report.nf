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
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_summary_*.tsv'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_site_bf_*.tsv'

    input:
    val  direction
    path json_files
    path fg_list_file  // optional: foreground species list (NO_FG_LIST sentinel when absent)

    output:
    path "FADE_report_${direction}.html",   emit: report
    path "fade_summary_${direction}.tsv",   emit: summary_tsv, optional: true
    path "fade_site_bf_${direction}.tsv",   emit: site_tsv,    optional: true

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

        n_json=\$(ls *.FADE.json 2>/dev/null | wc -l)
        echo "[FADE_REPORT] direction=${direction} | JSON files found: \${n_json}"
        if [ "\${n_json}" -eq 0 ]; then
            echo "[FADE_REPORT] WARNING: No *.FADE.json files found for direction '${direction}'. "\
                 "The report will contain no results. "\
                 "Possible causes: all FADE jobs failed, or no genes were selected for this direction."
        fi

        FADE_REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'FADE_report.Rmd',
                params = list(
                    json_dir        = '.',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}',
                    fg_list_file    = ${fg_arg}
                ),
                output_file = 'FADE_report_${direction}.html'
            )
        "

        if [ -f 'FADE_report_${direction}.html' ]; then
            echo "[FADE_REPORT] Report generated: FADE_report_${direction}.html"
        else
            echo "[FADE_REPORT] WARNING: Report file was not created for direction '${direction}'."
        fi
        """
    } else {
        """
        cp -R ${local_dir}/* .

        n_json=\$(ls *.FADE.json 2>/dev/null | wc -l)
        echo "[FADE_REPORT] direction=${direction} | JSON files found: \${n_json}"
        if [ "\${n_json}" -eq 0 ]; then
            echo "[FADE_REPORT] WARNING: No *.FADE.json files found for direction '${direction}'. "\
                 "The report will contain no results. "\
                 "Possible causes: all FADE jobs failed, or no genes were selected for this direction."
        fi

        FADE_REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                'FADE_report.Rmd',
                params = list(
                    json_dir        = '.',
                    direction       = '${direction}',
                    traitname       = '${traitname}',
                    bf_threshold    = ${bf_thr},
                    min_genes_hmap  = ${min_genes},
                    output_dir      = '${outdir}',
                    fg_list_file    = ${fg_arg}
                ),
                output_file = 'FADE_report_${direction}.html'
            )
        "

        if [ -f 'FADE_report_${direction}.html' ]; then
            echo "[FADE_REPORT] Report generated: FADE_report_${direction}.html"
        else
            echo "[FADE_REPORT] WARNING: Report file was not created for direction '${direction}'."
        fi
        """
    }
}
