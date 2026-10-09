#!/usr/bin/env nextflow
// ta_dataset_exploration.nf — Species-level exploration of the trait dataset (1.Dataset_exploration.Rmd).
// PhyloPhere | subworkflows/TRAIT_ANALYSIS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DATASET_EXPLORATION: renders 1.Dataset_exploration.Rmd, which plots the
 *  number of species per taxon (before and after pruning) and writes the
 *  per-species statistics table, trait_stats.csv (also read by
 *  2.Phenotype_exploration.Rmd). It runs in every configuration of the callers.
 *
 *  Consumes:  trait file, species tree, and the output directory of DATASET_PRUNE
 *             (or a NO_FILE sentinel when pruning was skipped)
 *  Produces:  data_exploration/ (figures, tables, trait_stats.csv, copies of the
 *             original trait and tree files) and the HTML report
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Dataset exploration report ───────────────────────────────────────────────

process DATASET_EXPLORATION {
    tag "dataset_exploration"
    label 'process_reporting_dataset'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path trait_file
    path tree_file
    path prune_results_dir

    output:
    path "data_exploration", emit: results_dir
    path "*.html", emit: reports, optional: true
    path "data_exploration/1.Data-exploration/1.Species_distribution/trait_stats.csv", emit: stats_file, optional: true
    path "data_exploration/**/*.csv", emit: data_tables, optional: true
    path "data_exploration/**/*.png", emit: plots, optional: true
    // Copies of the input trait and tree files. DATASET_PRUNE is optional, so
    // this is the one step that always leaves the original files in the results.
    path "data_exploration/1.Data-exploration/1.Species_distribution/original_trait_file.tsv", emit: original_trait_file, optional: true
    path "data_exploration/1.Data-exploration/1.Species_distribution/original_tree_file.nwk", emit: original_tree_file, optional: true

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
    def pss_top_pct = params.pss_top_pct ?: '0.05'
    def perm_strategy = params.perm_strategy ?: 'best_model'
    def trait_type = params.trait_type ?: ''
    def prune_dir = prune_results_dir ?: ''

    // The two branches are identical except that the container one runs Rscript through the image entrypoint.
    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        if [ -n "${prune_dir}" ] && [ -d "${prune_dir}" ] && [ "${prune_dir}" != "data_exploration" ]; then
          cp -R "${prune_dir}"/* data_exploration
        fi
        mkdir -p data_exploration/1.Data-exploration/1.Species_distribution
        cp '${trait_file}' data_exploration/1.Data-exploration/1.Species_distribution/original_trait_file.tsv
        cp '${tree_file}' data_exploration/1.Data-exploration/1.Species_distribution/original_tree_file.nwk
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '1.Dataset_exploration.Rmd',
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
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}'
                ),
                output_file = '1.Dataset_exploration.html',
                envir = new.env()
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        if [ -n "${prune_dir}" ] && [ -d "${prune_dir}" ] && [ "${prune_dir}" != "data_exploration" ]; then
          cp -R "${prune_dir}"/* data_exploration
        fi
        mkdir -p data_exploration/1.Data-exploration/1.Species_distribution
        cp '${trait_file}' data_exploration/1.Data-exploration/1.Species_distribution/original_trait_file.tsv
        cp '${tree_file}' data_exploration/1.Data-exploration/1.Species_distribution/original_tree_file.nwk
        Rscript -e "
            rmarkdown::render(
                '1.Dataset_exploration.Rmd',
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
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}'
                ),
                output_file = '1.Dataset_exploration.html',
                envir = new.env()
            )
        "
        """
    }
}
