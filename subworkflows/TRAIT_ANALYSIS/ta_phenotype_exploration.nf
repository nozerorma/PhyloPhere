#!/usr/bin/env nextflow

/*
#  Trait analysis: phenotype exploration (Rmarkdown)
*/

process PHENOTYPE_EXPLORATION {
    tag "phenotype_exploration"
    label 'process_reporting_phenotype'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path trait_file
    path tree_file
    path results_dir

    output:
    path "data_exploration"

    script:
    def local_dir = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local"
    def seed = params.seed ?: ''
    def clade = params.clade_name ?: ''
    def taxon = params.taxon_of_interest ?: ''
    def trait = params.traitname ?: ''
    def p_trait = params.p_trait ?: ''
    def n_trait = params.n_trait ?: ''
    def tax_id = params.tax_id ?: ''
    def branch_trait = params.branch_trait ?: ''
    def secondary_trait = params.secondary_trait ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        cp -R ${results_dir}/* data_exploration
        /usr/local/bin/_entrypoint.sh Rscript -e "rmarkdown::render('2.Phenotype_exploration.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          'data_exploration' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${p_trait}' \
          '${n_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}' 
        """
    } else {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        Rscript -e "rmarkdown::render('2.Phenotype_exploration.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          'data_exploration' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${p_trait}' \
          '${n_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}'
        """
    }
}
