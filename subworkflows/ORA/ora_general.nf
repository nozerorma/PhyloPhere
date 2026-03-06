#!/usr/bin/env nextflow

/*
#  ORA (WebGestalt) general process wrapper
*/

process ORA_GENERAL_REPORT {
    tag "ora_general"
    label 'process_reporting'

    publishDir path: "${params.outdir}/ora", mode: 'copy', overwrite: true, pattern: '{ORA_files/**,ora_results/**,ora_summary/**,ora_plots/**}'
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path background_file
    path gene_list_files

    output:
    path "*.html", emit: report
    path "ora_results/**", emit: ora_results, optional: true
    path "ora_summary/**", emit: ora_summary, optional: true
    path "ora_plots/**", emit: ora_plots, optional: true
    path "ORA_files/**", emit: assets, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ORA/local"
    def outdir = "${params.outdir}/ora"
    def project_name = params.ora_project_name ?: 'ORA_WebGestalt'
    def organism = params.ora_organism ?: 'hsapiens'
    def min_num = params.ora_min_num ?: 5
    def max_num = params.ora_max_num ?: 500
    def fdr_thr = params.ora_fdr ?: 0.1
    def report_num = params.ora_report_num ?: 20
    def top_thr = params.ora_top_thr ?: 15
    def n_threads = params.ora_threads ?: 8
    def enable_overlap = params.ora_enable_overlap_heatmap ? 'TRUE' : 'FALSE'
    def enable_enrichment_map = params.ora_enable_enrichment_map ? 'TRUE' : 'FALSE'
    def overlap_top_terms = params.ora_overlap_top_terms ?: 50
    def compare_metric = params.ora_compare_metric ?: 'overlap'
    def db_default = params.ora_databases_default ?: 'geneontology_Biological_Process_noRedundant,geneontology_Cellular_Component_noRedundant,geneontology_Molecular_Function_noRedundant,pathway_KEGG,pathway_Reactome,pathway_WikiPathways'
    def db_extra = params.ora_databases_extra ?: ''
    def db_combined = [db_default, db_extra].findAll { it && it.toString().trim() }.join(',')
    def bg_name = background_file.getName().replace("'", "\\'")
    // gene_list_files are staged as symlinks in the work directory; the Rmd
    // will auto-detect them via list.files(), so we no longer need to build
    // an explicit R vector here.  Keeping the variable to avoid unused-input
    // warnings from Nextflow.
    def gene_files_r = gene_list_files
        .collect { it.getName().replace("'", "\\'") }
        .collect { "'${it}'" }
        .join(', ')

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'ORA_general.Rmd',
                params = list(
                    background_file = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir = '${outdir}',
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
                output_file = 'ORA_general.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'ORA_general.Rmd',
                params = list(
                    background_file = '${background_file}',
                    background_basename = '${bg_name}',
                    output_dir = '${outdir}',
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
                output_file = 'ORA_general.html'
            )
        "
        """
    }
}
