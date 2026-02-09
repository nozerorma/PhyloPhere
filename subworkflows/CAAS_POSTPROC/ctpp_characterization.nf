#!/usr/bin/env nextflow

/*
#  CAAS Post-Processing: Characterization and reporting (Rmarkdown)
*/

process CAAS_POSTPROC_REPORT {
    tag "caas_postproc_report"
    label 'process_reporting'
    publishDir "${params.postproc_outdir}", mode: 'copy', overwrite: true
    
    input:
    path discovery_file
    path filter_summary
    path gene_ensembl_file
    
    output:
    path "reports", emit: report_dir
    path "reports/*.html", emit: reports, optional: true
    path "reports/*.tsv", emit: data_tables, optional: true
    path "reports/*.png", emit: plots, optional: true
    
    script:
    def local_dir = "${baseDir}/subworkflows/CAAS_POSTPROC/local"
    def discovery = discovery_file.toString()
    def filter_sum = filter_summary.toString()
    def gene_len = gene_ensembl_file.toString()
    def mode = params.caas_postproc_mode
    def outdir = params.postproc_outdir
    def extreme_thresh = params.extreme_threshold
    def iqr_mult = params.iqr_multiplier
    def gene_filter = params.gene_filter_mode
    def gen_manhattan = params.generate_manhattan
    def manhattan_min_dens = params.manhattan_min_density
    
    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        mkdir -p reports
        /usr/local/bin/_entrypoint.sh Rscript -e "rmarkdown::render('CAAS_postproc.Rmd', output_dir='reports', quiet=TRUE)" \
          '${discovery}' \
          '${filter_sum}' \
          '${gene_len}' \
          '${mode}' \
          '${outdir}' \
          '${extreme_thresh}' \
          '${iqr_mult}' \
          '${gene_filter}' \
          '${gen_manhattan}' \
          '${manhattan_min_dens}'
        """
    } else {
        """
        cp -R ${local_dir}/* .
        mkdir -p reports
        Rscript -e "rmarkdown::render('CAAS_postproc.Rmd', output_dir='reports', quiet=TRUE)" \
          '${discovery}' \
          '${filter_sum}' \
          '${gene_len}' \
          '${mode}' \
          '${outdir}' \
          '${extreme_thresh}' \
          '${iqr_mult}' \
          '${gene_filter}' \
          '${gen_manhattan}' \
          '${manhattan_min_dens}'
        """
    }
}
