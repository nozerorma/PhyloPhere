// CAAS Signification Processes
// Perform significance testing via hypergeometric and permutation analysis

process CAAS_SIGNIFICATION_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/signification", mode: 'copy', overwrite: true, pattern: '{CAAS_signification_files/**,gene_lists/**,meta_caas/**}'
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'
    
    input:
    path filtered_discovery
    path filtered_background
    path bootstrap_files
    
    output:
    path "*.html", emit: report
    path "CAAS_signification_files/**", emit: assets, optional: true
    path "gene_lists/**", emit: gene_lists, optional: true
    path "meta_caas/**", emit: meta_caas, optional: true
    
    script:
    def caap_mode_r = params.caap_mode ? 'TRUE' : 'FALSE'
    def local_dir = "${baseDir}/subworkflows/CAAS_SIGNIFICATION/local"
    // Bootstrap files are collected into a single list, get the first file (or all files)
    def bootstrap_file_name = bootstrap_files.getName()
    def outdir = "${params.outdir}/signification"


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'CAAS_signification.Rmd',
                params = list(
                    filtered_discovery = '${filtered_discovery}',
                    filtered_background = '${filtered_background}',
                    bootstrap_file = '${bootstrap_file_name}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = 'CAAS_signification.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        Rscript -e "
            rmarkdown::render(
                'CAAS_signification.Rmd',
                params = list(
                    filtered_discovery = '${filtered_discovery}',
                    filtered_background = '${filtered_background}',
                    bootstrap_file = '${bootstrap_file_name}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = 'CAAS_signification.html'
            )
        "
        """
    }
}