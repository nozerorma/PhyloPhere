#!/usr/bin/env nextflow

/*
#  Trait analysis: phenotype exploration (Rmarkdown)
*/

process PHENOTYPE_EXPLORATION {
    tag "phenotype_exploration"
    label 'process_reporting_phenotype'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path trait_file
    path tree_file
    path results_dir

    output:
    path "data_exploration", emit: results_dir
    path "*.html", emit: reports, optional: true
    path "data_exploration/**/*.csv", emit: data_tables, optional: true
    path "data_exploration/**/*.png", emit: plots, optional: true

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
    def discrete_method = params.discrete_method ?: 'quartile'
    def top_quantile = params.top_quantile ?: '0.75'
    def bottom_quantile = params.bottom_quantile ?: '0.25'
    def contrast_max_iter = params.contrast_max_iter ?: '3'

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        cp -R ${results_dir}/* data_exploration
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '2.Phenotype_exploration.Rmd',
                params = list(
                    trait_file = '${trait_file}',
                    tree_file = '${tree_file}',
                    output_dir = 'data_exploration',
                    seed = '${seed}',
                    clade_name = '${clade}',
                    taxon_of_interest = '${taxon}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    tax_id = '${tax_id}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}',
                    discrete_method = '${discrete_method}',
                    top_quantile = '${top_quantile}',
                    bottom_quantile = '${bottom_quantile}',
                    contrast_max_iter = '${contrast_max_iter}'
                ),
                output_file = '2.Phenotype_exploration.html',
                envir = new.env()
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "
            rmarkdown::render(
                '2.Phenotype_exploration.Rmd',
                params = list(
                    trait_file = '${trait_file}',
                    tree_file = '${tree_file}',
                    output_dir = 'data_exploration',
                    seed = '${seed}',
                    clade_name = '${clade}',
                    taxon_of_interest = '${taxon}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    tax_id = '${tax_id}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}',
                    discrete_method = '${discrete_method}',
                    top_quantile = '${top_quantile}',
                    bottom_quantile = '${bottom_quantile}',
                    contrast_max_iter = '${contrast_max_iter}'
                ),
                output_file = '2.Phenotype_exploration.html',
                envir = new.env()
            )
        "
        """
    }
}
