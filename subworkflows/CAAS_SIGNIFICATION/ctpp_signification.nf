// CAAS Signification Processes
// Perform significance testing via hypergeometric and permutation analysis

process CAAS_SIGNIFICATION_REPORT {
    label 'process_reporting'
    publishDir "${params.signification_outdir}", mode: 'copy', overwrite: true
    
    input:
    path filtered_discovery
    path filtered_background
    path bootstrap_files
    
    output:
    path "CAAS_signification.html", emit: report
    path "CAAS_signification_files/**", emit: assets, optional: true
    path "gene_lists/**", emit: gene_lists, optional: true
    path "commonalities/**", emit: commonalities, optional: true
    
    script:
    def caap_mode_r = params.caap_mode ? 'TRUE' : 'FALSE'
    def local_dir = "${baseDir}/subworkflows/CAAS_SIGNIFICATION/local"
    // Bootstrap files are collected into a single list, get the first file (or all files)
    def bootstrap_file_name = bootstrap_files.getName()
    """
    #!/usr/bin/env bash
    set -e
    
    # Copy Rmd file to working directory
    cp ${local_dir}/CAAS_signification.Rmd .
    
    # Render R Markdown report
    Rscript -e "
        rmarkdown::render(
            'CAAS_signification.Rmd',
            params = list(
                filtered_discovery = '${filtered_discovery}',
                filtered_background = '${filtered_background}',
                bootstrap_file = '${bootstrap_file_name}',
                output_dir = '.',
                caap_mode = ${caap_mode_r}
            ),
            output_file = 'CAAS_signification.html',
            envir = new.env()
        )
    "
    """
}
