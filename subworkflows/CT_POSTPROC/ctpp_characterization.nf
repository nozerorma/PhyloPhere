#!/usr/bin/env nextflow

/*
#  CT Post-Processing: Characterization and reporting (Rmarkdown)
*/

process CT_POSTPROC_REPORT {
    tag "caas_postproc_report"
    label 'process_reporting'
    publishDir path: "${params.outdir}/postproc", mode: 'copy', overwrite: true, pattern: '{CT_postproc_files/**,outliers/**,clusters/**,summary_statistics/**,disambiguation_characterization/**}'
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'
    
    input:
    path discovery_file
    path filter_summary
    val filter_dir
    path gene_ensembl_file
    path gene_stats_file
    
    output:
    path "*.html", emit: report
    path "CT_postproc_files/**", emit: assets, optional: true
    path "outliers/**", emit: outliers, optional: true
    path "clusters/**", emit: clusters, optional: true
    path "summary_statistics/**", emit: summary_stats, optional: true
    path "disambiguation_characterization/**", emit: disambiguation_characterization, optional: true
    
    script:
    def local_dir = "${baseDir}/subworkflows/CT_POSTPROC/local"
    def discovery = discovery_file.toString()
    def filter_sum = filter_summary.toString()
    def filter_dir_path = filter_dir
    def gene_len = gene_ensembl_file.toString()
    def gene_stats = gene_stats_file.toString()
    def mode = params.caas_postproc_mode
    def outdir = "${params.outdir}/postproc"
    def extreme_thresh = params.extreme_threshold
    def iqr_mult = params.iqr_multiplier
    def alpha_thresh = params.alpha_threshold
    def gene_filter = params.gene_filter_mode


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true
        
        REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'CT_postproc.Rmd',
                params = list(
                    discovery_file = '${discovery}',
                    filter_summary_file = '${filter_sum}',
                    filter_dir = '${filter_dir_path}',
                    gene_ensembl_file = '${gene_len}',
                    gene_stats_file = '${gene_stats}',
                    processing_mode = '${mode}',
                    output_dir = '${outdir}',
                    extreme_threshold = ${extreme_thresh},
                    iqr_multiplier = ${iqr_mult},
                    alpha_threshold = ${alpha_thresh},
                    gene_filter_mode = '${gene_filter}'
                ),
                output_file = 'CT_postproc.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true
        
        REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                'CT_postproc.Rmd',
                params = list(
                    discovery_file = '${discovery}',
                    filter_summary_file = '${filter_sum}',
                    filter_dir = '${filter_dir_path}',
                    gene_ensembl_file = '${gene_len}',
                    gene_stats_file = '${gene_stats}',
                    processing_mode = '${mode}',
                    output_dir = '${outdir}',
                    extreme_threshold = ${extreme_thresh},
                    iqr_multiplier = ${iqr_mult},
                    alpha_threshold = ${alpha_thresh},
                    gene_filter_mode = '${gene_filter}'
                ),
                output_file = 'CT_postproc.html'
            )
        "
        """
    }
}
