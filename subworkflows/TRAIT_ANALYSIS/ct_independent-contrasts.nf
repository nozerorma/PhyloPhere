#!/usr/bin/env nextflow

/*
#  Contrast selection for CT analysis (Rmarkdown)
*/

process CONTRAST_ALGORITHM {
    tag "CONTRAST_ALGORITHM"
    label 'process_contrast_selection'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path trait_file
    path tree_file
    path results_dir

    output:
    path "data_exploration", emit: contrast_results_dir
    path "data_exploration/2.CAAS/1.Traitfiles/traitfile.tab", emit: trait_file_out
    path "data_exploration/2.CAAS/3.Tree/pruned_tree_file.nwk", emit: tree_file_out
    path "data_exploration/2.CAAS/2.Bootstrap_traitfiles/boot_traitfile.tab", emit: bootstrap_trait_file_out

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
        mkdir -p ${results_dir}/HTML_reports
        /usr/local/bin/_entrypoint.sh Rscript -e "rmarkdown::render('4.Independent_contrasts.Rmd', output_dir='${results_dir}/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          '${results_dir}' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${n_trait}' \
          '${c_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}'
        
        # Copy the tree file to the expected output location
        mkdir -p ${results_dir}/2.CAAS/3.Tree
        cp ${tree_file} ${results_dir}/2.CAAS/3.Tree/pruned_tree_file.nwk
        """
    } else {
        """
        cp -R ${local_dir}/* .
        mkdir -p ${results_dir}/HTML_reports
        Rscript -e "rmarkdown::render('4.Independent_contrasts.Rmd', output_dir='${results_dir}/HTML_reports', quiet=TRUE)" --args \
          '${trait_file}' \
          '${tree_file}' \
          '${results_dir}' \
          '${seed}' \
          '${clade}' \
          '${taxon}' \
          '${trait}' \
          '${n_trait}' \
          '${c_trait}' \
          '${tax_id}' \
          '${secondary_trait}' \
          '${branch_trait}'
        
        # Copy the tree file to the expected output location
        mkdir -p ${results_dir}/2.CAAS/3.Tree
        cp ${tree_file} ${results_dir}/2.CAAS/3.Tree/pruned_tree_file.nwk
        """
    }
}
