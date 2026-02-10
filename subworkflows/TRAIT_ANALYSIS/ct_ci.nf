#!/usr/bin/env nextflow

/*
#  Contrast selection for CT analysis (Rmarkdown)
*/

process CI_COMPOSITION_REPORT {
    tag "CI_COMPOSITION_REPORT"
    label 'process_contrast_selection'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path trait_file
    path tree_file
    path results_dir

    output:
    path results_dir, emit: results_dir
    path "*.html", emit: reports, optional: true
    path "${results_dir}/**/*.csv", emit: data_tables, optional: true
    path "${results_dir}/**/*.png", emit: plots, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local"
    def seed = params.seed ?: ''
    def clade = params.clade_name ?: ''
    def taxon = params.taxon_of_interest ?: ''
    def trait = params.traitname ?: ''
    def n_trait = params.n_trait ?: ''
    def c_trait = params.c_trait ?: ''
    def tax_id = params.tax_id ?: ''
    def branch_trait = params.branch_trait ?: ''
    def secondary_trait = params.secondary_trait ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '3.CI-composition.Rmd',
                params = list(
                    trait_file = '${trait_file}',
                    tree_file = '${tree_file}',
                    output_dir = '${results_dir}',
                    seed = '${seed}',
                    clade_name = '${clade}',
                    taxon_of_interest = '${taxon}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    tax_id = '${tax_id}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}'
                ),
                output_file = '3.CI-composition.html',
                envir = new.env()
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "
            rmarkdown::render(
                '3.CI-composition.Rmd',
                params = list(
                    trait_file = '${trait_file}',
                    tree_file = '${tree_file}',
                    output_dir = '${results_dir}',
                    seed = '${seed}',
                    clade_name = '${clade}',
                    taxon_of_interest = '${taxon}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    tax_id = '${tax_id}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}'
                ),
                output_file = '3.CI-composition.html',
                envir = new.env()
            )
        "
        """
    }
}
