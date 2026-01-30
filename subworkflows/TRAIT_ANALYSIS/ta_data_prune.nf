#!/usr/bin/env nextflow

/*
#  Trait analysis: optional data pruning (Rmarkdown)
*/

process DATASET_PRUNE {
    tag "dataset_prune"
    label 'process_reporting_dataset'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path trait_file
    path tree_file

    output:
    path "data_exploration", emit: pruned_results_dir
    path "data_exploration/0.Data-pruning/pruned_trait_file.tsv", emit: pruned_trait_file
    path "data_exploration/0.Data-pruning/pruned_tree_file.nwk", emit: pruned_tree_file
    path "data_exploration/0.Data-pruning/pruned_trait_stats.csv", emit: pruned_stats_file

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
    def prune_list = params.prune_list ?: ''
    def prune_list_secondary = params.prune_list_secondary ?: ''


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        /usr/local/bin/_entrypoint.sh Rscript -e "rmarkdown::render('0.Data_pruning.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
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
          '${branch_trait}' \
          '${prune_list}' \
          '${prune_list_secondary}'

        """
    } else {
        """
        cp -R ${local_dir}/* .
        mkdir -p data_exploration/HTML_reports
        Rscript -e "rmarkdown::render('0.Data_pruning.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args \
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
          '${branch_trait}' \
          '${prune_list}' \
          '${prune_list_secondary}'
        """
    }
}
