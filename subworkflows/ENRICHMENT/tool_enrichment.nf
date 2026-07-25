#!/usr/bin/env nextflow

/*
 * TOOL_STRING_REPORT / TOOL_STRING_MODULE_REPORT
 * ───────────────────────────────────────────────
 * Per-tool STRING enrichment. Dynamic publishDir driven by `output_subpath`
 * so results land in the tool's own directory (e.g. RERConverge/string/,
 * FADE/string/, CT_ACCUMULATION/string/).
 *
 * TOOL_STRING_MODULE_REPORT wraps DOMINO_MODULES (subworkflows/ENRICHMENT/domino.nf)
 * ahead of TOOL_STRING_REPORT under the exact same call signature, so the three
 * existing call sites (RER_STRING_REPORT / FADE_STRING_REPORT /
 * ACCUMULATION_STRING_REPORT in rerconverge.nf / fade.nf / ct_accumulation.nf)
 * only need their `include` line repointed at this workflow, not their calls.
 *
 * Inputs (both)
 * ──────
 *   output_subpath  : val  — subpath under params.outdir
 *   background_file : path — all-genes-tested background .txt
 *   gene_list_files : path — one or more gene list .txt files
 *   report_label    : val  — filename stem for the output HTML
 */

include { DOMINO_MODULES } from './domino.nf'

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
    path domino_network_sif
    path domino_modules_dir

    output:
    path "*.html",              emit: report
    path "string_results/**",   emit: string_results,  optional: true
    path "string_summary/**",   emit: string_summary,  optional: true
    path "string_plots/**",     emit: string_plots,    optional: true
    path "string_networks/**",  emit: string_networks, optional: true

    script:
    def local_dir         = "${baseDir}/subworkflows/ENRICHMENT/local"
    def outdir            = "${params.outdir}/${output_subpath}/string"
    def project_name      = "STRING_${params.traitname ?: 'unknown_trait'}"
    def species           = params.string_species           ?: 9606
    def required_score    = params.string_required_score    ?: 400
    def net_score         = params.domino_network_score_thr ?: 700
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
                '15.STRING_report.Rmd',
                params = list(
                    background_file    = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir         = '${outdir}',
                    project_name       = '${report_label}',
                    species            = ${species},
                    required_score     = ${required_score},
                    domino_network_score_thr = ${net_score},
                    fdr_thr            = ${fdr_thr},
                    top_thr            = ${top_thr},
                    report_num         = ${report_num},
                    string_db_dir      = '${params.string_db_dir}',
                    domino_network_sif = '${domino_network_sif}',
                    domino_modules_dir = '${domino_modules_dir}'
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
                '15.STRING_report.Rmd',
                params = list(
                    background_file    = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir         = '${outdir}',
                    project_name       = '${report_label}',
                    species            = ${species},
                    required_score     = ${required_score},
                    domino_network_score_thr = ${net_score},
                    fdr_thr            = ${fdr_thr},
                    top_thr            = ${top_thr},
                    report_num         = ${report_num},
                    string_db_dir      = '${params.string_db_dir}',
                    domino_network_sif = '${domino_network_sif}',
                    domino_modules_dir = '${domino_modules_dir}'
                ),
                output_file = '${report_label}.html'
            )
        "
        """
    }
}

workflow TOOL_STRING_MODULE_REPORT {
    take:
    output_subpath
    background_file
    gene_list_files
    report_label

    main:
    def domino_out = DOMINO_MODULES(
        background_file,
        gene_list_files,
        params.domino_network_score_thr ?: 700,
        params.domino_slice_thr ?: 0.3,
        params.domino_module_thr ?: 0.05,
        params.string_db_dir ?: ''
    )
    TOOL_STRING_REPORT(
        output_subpath, background_file, gene_list_files, report_label,
        domino_out.network_sif, domino_out.modules_dir
    )

    emit:
    report          = TOOL_STRING_REPORT.out.report
    string_results  = TOOL_STRING_REPORT.out.string_results
    string_summary  = TOOL_STRING_REPORT.out.string_summary
    string_plots    = TOOL_STRING_REPORT.out.string_plots
    string_networks = TOOL_STRING_REPORT.out.string_networks
}
