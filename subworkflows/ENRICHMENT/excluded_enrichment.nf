#!/usr/bin/env nextflow

/*
#  Enrichment (clusterProfiler) excluded-genes process wrapper
#  Mirrors enrichment_general.nf but publishes to ${params.outdir}/enrichment_excluded
#  so that excluded-gene enrichment results do not overwrite standard enrichment outputs.
*/

process EXCLUDED_ENRICHMENT_REPORT {
    tag "enrichment_excluded_general"
    label 'process_reporting'

    publishDir path: "${params.outdir}/enrichment_excluded", mode: 'copy', overwrite: true, pattern: '{ORA_files/**,ora_results/**,ora_summary/**,ora_plots/**}'
    publishDir path: "${params.outdir}/html_reports",  mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path background_file
    path gene_list_files

    output:
    path "*.html",         emit: report
    path "ora_results/**", emit: ora_results, optional: true
    path "ora_summary/**", emit: ora_summary, optional: true
    path "ora_plots/**",   emit: ora_plots,   optional: true
    path "ORA_files/**",   emit: assets,      optional: true

    script:
    def local_dir    = "${baseDir}/subworkflows/ENRICHMENT/local"
    def outdir       = "${params.outdir}/enrichment_excluded"
    def gmt_dir      = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def project_name = "Enrichment_Excluded_${params.traitname ?: 'unknown_trait'}"
    def organism     = params.enrichment_organism     ?: 'hsapiens'
    def min_num      = params.enrichment_min_num      ?: 5
    def max_num      = params.enrichment_max_num      ?: 500
    def fdr_thr      = params.enrichment_fdr          ?: 0.1
    def report_num   = params.enrichment_report_num   ?: 20
    def top_thr      = params.enrichment_top_thr      ?: 15
    def n_threads    = task.cpus ?: 8
    def enable_overlap        = 'TRUE'
    def enable_enrichment_map = 'TRUE'
    def overlap_top_terms     = params.enrichment_overlap_top_terms ?: 50
    def compare_metric        = params.enrichment_compare_metric    ?: 'overlap'
    // ora_databases is vestigial: 14.ORA_report.Rmd loads gene sets from gmt_dir and only
    // uses this list for a non-empty guard + a cosmetic printout. Kept as a fixed label.
    def db_combined = 'gmt'
    def bg_name     = background_file.getName().replace("'", "\\'")
    def gene_files_r = gene_list_files
        .collect { it.getName().replace("'", "\\'") }
        .collect { "'${it}'" }
        .join(', ')

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '14.ORA_report.Rmd',
                params = list(
                    background_file = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir = '${outdir}',
                    gmt_dir = '${gmt_dir}',
                    project_name = '${project_name}',
                    organism = '${organism}',
                    ora_databases = '${db_combined}',
                    min_num = ${min_num},
                    max_num = ${max_num},
                    fdr_thr = ${fdr_thr},
                    report_num = ${report_num},
                    top_thr = ${top_thr},
                    n_threads = ${n_threads},
                    enable_overlap_heatmap = ${enable_overlap},
                    enable_enrichment_map = ${enable_enrichment_map},
                    overlap_top_terms = ${overlap_top_terms},
                    compare_metric = '${compare_metric}'
                ),
                output_file = '14.Enrichment_excluded.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                '14.ORA_report.Rmd',
                params = list(
                    background_file = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir = '${outdir}',
                    gmt_dir = '${gmt_dir}',
                    project_name = '${project_name}',
                    organism = '${organism}',
                    ora_databases = '${db_combined}',
                    min_num = ${min_num},
                    max_num = ${max_num},
                    fdr_thr = ${fdr_thr},
                    report_num = ${report_num},
                    top_thr = ${top_thr},
                    n_threads = ${n_threads},
                    enable_overlap_heatmap = ${enable_overlap},
                    enable_enrichment_map = ${enable_enrichment_map},
                    overlap_top_terms = ${overlap_top_terms},
                    compare_metric = '${compare_metric}'
                ),
                output_file = '14.Enrichment_excluded.html'
            )
        "
        """
    }
}
