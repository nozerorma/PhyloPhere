#!/usr/bin/env nextflow

/*
#  Trait analysis: optional data pruning (Rmarkdown)
*/

process DATASET_PRUNE {
    tag "dataset_prune"
    label 'process_reporting_dataset'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/HTML_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path trait_file
    path tree_file

    output:
    path "data_exploration", emit: pruned_results_dir
    path "*.html", emit: reports, optional: true
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
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '0.Data_pruning.Rmd',
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
                    prune_list = '${prune_list}',
                    prune_list_secondary = '${prune_list_secondary}'
                ),
                output_file = '0.Data_pruning.html',
                envir = new.env()
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "
            rmarkdown::render(
                '0.Data_pruning.Rmd',
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
                    prune_list = '${prune_list}',
                    prune_list_secondary = '${prune_list_secondary}'
                ),
                output_file = '0.Data_pruning.html',
                envir = new.env()
            )
        "
        """
    }
}
