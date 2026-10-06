#!/usr/bin/env nextflow
// ct_ci.nf — Render the trait composition report that prepares the candidate contrast pairs.
// PhyloPhere | subworkflows/TRAIT_ANALYSIS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CI_COMPOSITION_REPORT: renders 3.CI-composition.Rmd. Depending on the trait it
 *  builds Jeffreys credible intervals (count columns), uses the coded levels (ordinal)
 *  or scores species pairs by Phylogenetic Shift Score (continuous), and writes the
 *  pairwise table and the candidate foreground/background species pool.
 *
 *  Consumes:  trait file, species tree, results directory of the previous exploration step
 *  Produces:  the updated results directory, html_reports/*.html, candidate_species.tab
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Composition report ─────────────────────────────────────────────────────────

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
    // Candidate foreground/background species pool before the Dunn-based selection
    // (candidate_species.tab format), read by SELECTION_PREP for FADE.
    path "${results_dir}/1.Data-exploration/5.CI_overlaps/candidate_species.tab", emit: candidate_species_out, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local"
    // The shared contrast-selection core (lean_contrast_selector.R) and the PSS engine
    // (pss_core.R) live with the CT scripts; they are copied into src/ so the Rmd and
    // selection_algorithm.R find them through the getwd()/src/ path.
    def ct_scripts = "${baseDir}/subworkflows/CT/local/scripts"
    def seed = params.seed ?: ''
    def clade = params.clade_name ?: ''
    def taxon = params.taxon_of_interest ?: ''
    def sp_colname = params.sp_colname ?: 'species'
    def trait = params.traitname ?: ''
    def n_trait = params.n_trait ?: ''
    def c_trait = params.c_trait ?: ''
    def tax_id = params.tax_id ?: ''
    def branch_trait = params.branch_trait ?: ''
    def secondary_trait = params.secondary_trait ?: ''
    def pss_top_pct = params.pss_top_pct ?: '0.05'
    def perm_strategy = params.perm_strategy ?: 'best_model'
    def trait_type = params.trait_type ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        cp ${ct_scripts}/lean_contrast_selector.R ${ct_scripts}/pss_core.R src/
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
                    sp_colname = '${sp_colname}',
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
        cp ${ct_scripts}/lean_contrast_selector.R ${ct_scripts}/pss_core.R src/
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
                    sp_colname = '${sp_colname}',
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
