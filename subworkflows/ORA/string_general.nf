#!/usr/bin/env nextflow

/*
#  STRING (rbioapi) enrichment general process wrapper
*/

process STRING_GENERAL_REPORT {
    tag "string_general"
    label 'process_reporting'

    publishDir path: "${params.outdir}/string", mode: 'copy', overwrite: true, pattern: '{string_results/**,string_summary/**,string_plots/**}'
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path background_file
    path gene_list_files

    output:
    path "*.html",            emit: report
    path "string_results/**", emit: string_results, optional: true
    path "string_summary/**", emit: string_summary, optional: true
    path "string_plots/**",   emit: string_plots,   optional: true

    script:
    def local_dir              = "${baseDir}/subworkflows/ORA/local"
    def outdir                 = "${params.outdir}/string"
    def project_name           = params.string_project_name      ?: 'STRING_Analysis'
    def species                = params.string_species            ?: 9606
    def required_score         = params.string_required_score     ?: 400
    def fdr_thr                = params.string_fdr               ?: 0.1
    def top_thr                = params.string_top_thr           ?: 15
    def report_num             = params.string_report_num        ?: 20
    def enable_ppi             = params.string_enable_ppi        ? 'TRUE' : 'FALSE'
    def enable_overlap         = params.string_enable_overlap_heatmap  ? 'TRUE' : 'FALSE'
    def enable_enrichment_map  = params.string_enable_enrichment_map   ? 'TRUE' : 'FALSE'
    def overlap_top_terms      = params.string_overlap_top_terms ?: 50
    def compare_metric         = params.string_compare_metric    ?: 'overlap'
    def bg_name                = background_file.getName().replace("'", "\\'")
    def gene_files_r           = gene_list_files
        .collect { it.getName().replace("'", "\\'") }
        .collect { "'${it}'" }
        .join(', ')

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'STRING_general.Rmd',
                params = list(
                    background_file       = '${background_file}',
                    gene_list_files       = c(${gene_files_r}),
                    background_basename   = '${bg_name}',
                    output_dir            = '${outdir}',
                    project_name          = '${project_name}',
                    species               = ${species},
                    required_score        = ${required_score},
                    fdr_thr               = ${fdr_thr},
                    top_thr               = ${top_thr},
                    report_num            = ${report_num},
                    enable_ppi_enrichment = ${enable_ppi},
                    enable_overlap_heatmap   = ${enable_overlap},
                    enable_enrichment_map    = ${enable_enrichment_map},
                    overlap_top_terms     = ${overlap_top_terms},
                    compare_metric        = '${compare_metric}'
                ),
                output_file = 'STRING_general.html'
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
                    background_file       = '${background_file}',
                    gene_list_files       = c(${gene_files_r}),
                    background_basename   = '${bg_name}',
                    output_dir            = '${outdir}',
                    project_name          = '${project_name}',
                    species               = ${species},
                    required_score        = ${required_score},
                    fdr_thr               = ${fdr_thr},
                    top_thr               = ${top_thr},
                    report_num            = ${report_num},
                    enable_ppi_enrichment = ${enable_ppi},
                    enable_overlap_heatmap   = ${enable_overlap},
                    enable_enrichment_map    = ${enable_enrichment_map},
                    overlap_top_terms     = ${overlap_top_terms},
                    compare_metric        = '${compare_metric}'
                ),
                output_file = 'STRING_general.html'
            )
        "
        """
    }
}
