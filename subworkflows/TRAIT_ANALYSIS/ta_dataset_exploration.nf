#!/usr/bin/env nextflow

/*
#  Trait analysis: dataset exploration (Rmarkdown)
*/

process DATASET_EXPLORATION {
    tag "dataset_exploration"
    label 'process_reporting_dataset'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path trait_file
    path tree_file
    path prune_results_dir

    output:
    path "data_exploration"

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
    def prune_dir = prune_results_dir ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        if [ -n "${prune_dir}" ] && [ -d "${prune_dir}" ] && [ "${prune_dir}" != "data_exploration" ]; then
          cp -R "${prune_dir}"/* data_exploration
        fi
        /usr/local/bin/_entrypoint.sh Rscript -e "rmarkdown::render('1.Dataset_exploration.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          'data_exploration' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${n_trait}'  \
          '${c_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}'

        """
    } else {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        if [ -n "${prune_dir}" ] && [ -d "${prune_dir}" ] && [ "${prune_dir}" != "data_exploration" ]; then
          cp -R "${prune_dir}"/* data_exploration
        fi
        Rscript -e "rmarkdown::render('1.Dataset_exploration.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          'data_exploration' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${n_trait}' \
          '${c_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}'
        """
    }
}
