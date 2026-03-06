// CT Signification Processes
// Perform significance testing via hypergeometric and permutation analysis

process CAAS_SIGNIFICATION_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/signification", mode: 'copy', overwrite: true, pattern: '{CAAS_signification_files/**,gene_lists/**,meta_caas/**}'
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'
    
    input:
    path discovery_input
    path background_input
    path bootstrap_files
    
    output:
    path "*.html", emit: report
    path "CAAS_signification_files/**", emit: assets, optional: true
    path "gene_lists/**", emit: gene_lists, optional: true
    path "meta_caas/**", emit: meta_caas, optional: true
    path "meta_caas/global_meta_caas.tsv", emit: global_meta_caas, optional: true
    
    script:
    def caap_mode_r = params.caap_mode ? 'TRUE' : 'FALSE'
    def local_dir = "${baseDir}/subworkflows/CT_SIGNIFICATION/local"
    // Bootstrap files are collected into a single list, get the first file (or all files)
    def bootstrap_input = bootstrap_files.getName()
    def outdir = "${params.outdir}/signification"


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'CT_signification.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    bootstrap_file = '${bootstrap_input}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = 'CT_signification.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        Rscript -e "
            rmarkdown::render(
                'CT_signification.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    bootstrap_file = '${bootstrap_input}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = 'CT_signification.html'
            )
        "
        """
    }
}
