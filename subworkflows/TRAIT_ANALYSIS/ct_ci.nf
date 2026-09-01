#!/usr/bin/env nextflow

/*
#  Contrast selection for CT analysis (Rmarkdown)
*/

process CI_COMPOSITION_REPORT {
    tag "CI_COMPOSITION_REPORT"
    label 'process_contrast_selection'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.equals('data_exploration') || filename.startsWith('data_exploration/') ? filename : null }
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path trait_file
    path tree_file
    path results_dir

    output:
    path results_dir, emit: results_dir
    path "*.html", emit: reports, optional: true
    path "${results_dir}/**/*.csv", emit: data_tables, optional: true
    path "${results_dir}/**/*.png", emit: plots, optional: true
    // Pre-Dunn candidate foreground/background species pool (traitfile format),
    // consumed by SELECTION_PREP -> EXTRACT_EXTREME_SPECIES for FADE.
    path "${results_dir}/1.Data-exploration/5.CI_overlaps/candidate_species.tab", emit: candidate_species_out, optional: true

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
    def pss_top_pct = params.pss_top_pct ?: '0.01'
    def perm_strategy = params.perm_strategy ?: 'best_model'
    def trait_type = params.trait_type ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        if [ -L "${results_dir}" ]; then
            target=\$(readlink -f "${results_dir}")
            rm -f "${results_dir}"
            cp -r "\${target}" "${results_dir}"
        fi
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
                    branch_trait = '${branch_trait}',
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}'
                ),
                output_file = '3.CI-composition.html',
                envir = new.env()
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        if [ -L "${results_dir}" ]; then
            target=\$(readlink -f "${results_dir}")
            rm -f "${results_dir}"
            cp -r "\${target}" "${results_dir}"
        fi
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
                    branch_trait = '${branch_trait}',
                    trait_type = '${trait_type}',
                    pss_top_pct = '${pss_top_pct}',
                    perm_strategy = '${perm_strategy}'
                ),
                output_file = '3.CI-composition.html',
                envir = new.env()
            )
        "
        """
    }
}
