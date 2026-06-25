#!/usr/bin/env nextflow

/*
 * TOOL_STRING_REPORT
 * ──────────────────
 * Per-tool STRING enrichment. Dynamic publishDir driven by `output_subpath`
 * so results land in the tool's own directory (e.g. RERConverge/string/,
 * FADE/string/, CT_ACCUMULATION/string/).
 *
 * Inputs
 * ──────
 *   output_subpath  : val  — subpath under params.outdir
 *   background_file : path — all-genes-tested background .txt
 *   gene_list_files : path — one or more gene list .txt files
 *   report_label    : val  — filename stem for the output HTML
 */

process TOOL_STRING_REPORT {
    tag "tool_string|${report_label}"
    label 'process_reporting'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    publishDir path: { "${params.outdir}/${output_subpath.toLowerCase()}/string" },
               mode: 'copy', overwrite: true,
               pattern: '{string_results/**,string_summary/**,string_plots/**,string_networks/**}'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    val  output_subpath
    path background_file
    path gene_list_files
    val  report_label

    output:
    path "*.html",              emit: report
    path "string_results/**",   emit: string_results,  optional: true
    path "string_summary/**",   emit: string_summary,  optional: true
    path "string_plots/**",     emit: string_plots,    optional: true
    path "string_networks/**",  emit: string_networks, optional: true

    script:
    def local_dir         = "${baseDir}/subworkflows/ENRICHMENT/local"
    def outdir            = "${params.outdir}/${output_subpath}/string"
    def project_name      = params.string_project_name     ?: 'STRING_Analysis'
    def species           = params.string_species           ?: 9606
    def required_score    = params.string_required_score    ?: 400
    def network_score_thr = params.string_network_score_thr ?: 700
    def mcl_inflation     = params.string_mcl_inflation     ?: 2.5
    def fdr_thr           = params.string_fdr              ?: 0.1
    def top_thr           = params.string_top_thr          ?: 15
    def report_num        = params.string_report_num       ?: 20
    def bg_name           = background_file instanceof List
        ? background_file[0].getName().replace("'", "\\'")
        : background_file.getName().replace("'", "\\'")

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'STRING_general.Rmd',
                params = list(
                    background_file    = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir         = '${outdir}',
                    project_name       = '${report_label}',
                    species            = ${species},
                    required_score     = ${required_score},
                    network_score_thr  = ${network_score_thr},
                    mcl_inflation      = ${mcl_inflation},
                    fdr_thr            = ${fdr_thr},
                    top_thr            = ${top_thr},
                    report_num         = ${report_num}
                ),
                output_file = '${report_label}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'STRING_general.Rmd',
                params = list(
                    background_file    = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir         = '${outdir}',
                    project_name       = '${report_label}',
                    species            = ${species},
                    required_score     = ${required_score},
                    network_score_thr  = ${network_score_thr},
                    mcl_inflation      = ${mcl_inflation},
                    fdr_thr            = ${fdr_thr},
                    top_thr            = ${top_thr},
                    report_num         = ${report_num}
                ),
                output_file = '${report_label}.html'
            )
        "
        """
    }
}
