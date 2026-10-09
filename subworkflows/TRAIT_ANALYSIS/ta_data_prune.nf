#!/usr/bin/env nextflow
// ta_data_prune.nf — Remove listed and phenotype-less species from the trait table and tree.
// PhyloPhere | subworkflows/TRAIT_ANALYSIS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DATASET_PRUNE: renders 0.Data_pruning.Rmd, which drops the species of
 *  params.prune_list (required) and the species with a missing phenotype or
 *  count value from the trait table, the tree and the per-species statistics.
 *  Species of params.prune_list_secondary only lose their secondary-trait value
 *  in the plots. Runs only when params.prune_data is set.
 *
 *  Consumes:  trait file, species tree
 *  Produces:  data_exploration/0.Data-pruning/ (pruned_trait_file.tsv,
 *             pruned_tree_file.nwk, pruned_trait_stats.csv), the report figures
 *             and tables under data_exploration/, and the HTML report
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Pruning report ───────────────────────────────────────────────────────────

process DATASET_PRUNE {
    tag "dataset_prune"
    label 'process_reporting_dataset'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

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
    def sp_colname = params.sp_colname ?: 'species'
    def trait = params.traitname ?: ''
    def n_trait = params.n_trait ?: ''
    def c_trait = params.c_trait ?: ''
    def branch_trait = params.branch_trait ?: ''
    def secondary_trait = params.secondary_trait ?: ''
    def prune_list = params.prune_list ?: ''
    def prune_list_secondary = params.prune_list_secondary ?: ''
    def pss_top_pct = params.pss_top_pct ?: '0.05'
    def perm_strategy = params.perm_strategy ?: 'best_model'
    def trait_type = params.trait_type ?: ''
    def max_contrasts = params.max_contrasts ?: '0'

    // The two branches are identical except that the container one runs Rscript through the image entrypoint.
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
                    sp_colname = '${sp_colname}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}',
                    prune_list = '${prune_list}',
                    prune_list_secondary = '${prune_list_secondary}',
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}',
                    max_contrasts = '${max_contrasts}'
                ),
                output_file = '2.Phenotype_exploration_pruned.html',
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
                    sp_colname = '${sp_colname}',
                    traitname = '${trait}',
                    n_trait = '${n_trait}',
                    c_trait = '${c_trait}',
                    secondary_trait = '${secondary_trait}',
                    branch_trait = '${branch_trait}',
                    prune_list = '${prune_list}',
                    prune_list_secondary = '${prune_list_secondary}',
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}',
                    max_contrasts = '${max_contrasts}'
                ),
                output_file = '2.Phenotype_exploration_pruned.html',
                envir = new.env()
            )
        "
        """
    }
}
