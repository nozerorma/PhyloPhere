#!/usr/bin/env nextflow

/*
 * FADE_RUN
 * ────────
 * Run HyPhy FADE (FUBAR Approach to Directional Evolution) on a single gene
 * alignment + annotated tree, testing foreground branches for directional
 * amino-acid selection.
 *
 * Inputs
 * ──────
 *   gene_id        : string — gene identifier (used for file naming/tagging)
 *   direction      : string — 'top' or 'bottom'
 *   fasta          : protein FASTA alignment (taxa names must match tree labels)
 *   annotated_tree : Newick tree with {Foreground} leaf labels (from ANNOTATE_TREE_FG)
 *
 * Outputs
 * ───────
 *   fade_json : FADE results JSON (HyPhy standard output)
 */

process FADE_RUN {
    tag "${gene_id}|${direction}"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/fade/${direction}/json",
               mode: 'copy', overwrite: true,
               pattern: '*.FADE.json'

    errorStrategy 'ignore'  // Skip genes that fail (e.g., too few FG branches)

    input:
    tuple val(gene_id), val(direction), path(fasta), path(annotated_tree)
    path lg_dat

    output:
    tuple val(gene_id), val(direction), path("${gene_id}.FADE.json"), emit: fade_json, optional: true

    script:
    def model   = params.fade_model   ?: 'LG'
    def model_file_arg = model == 'LG' ? "--model-file lg.dat" : ""
    def method  = params.fade_method  ?: 'Variational-Bayes'
    def grid    = params.fade_grid    ?: 20
    def conc    = params.fade_concentration ?: 0.5

    // MCMC-specific args (only meaningful when method != Variational-Bayes)
    def mcmc_args = (method == 'Variational-Bayes') ? "" :
        """--chains ${params.fade_chains ?: 5} \\
           --chain-length ${params.fade_chain_length ?: 2000000} \\
           --burn-in ${params.fade_burn_in ?: 1000000} \\
           --samples ${params.fade_samples ?: 1000}"""

    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh hyphy fade \\
            --alignment "${fasta}" \\
            --tree      "${annotated_tree}" \\
            --branches  Foreground \\
            --model     ${model} \\
            ${model_file_arg} \\
            --method    "${method}" \\
            --grid      ${grid} \\
            --concentration_parameter ${conc} \\
            ${mcmc_args} \\
            --output    "${gene_id}.FADE.json" \\
        || echo "FADE failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.FADE.json" ] || rm -f "${gene_id}.FADE.json"
        """
    } else {
        """
        hyphy fade \\
            --alignment "${fasta}" \\
            --tree      "${annotated_tree}" \\
            --branches  Foreground \\
            --model     ${model} \\
            ${model_file_arg} \\
            --method    "${method}" \\
            --grid      ${grid} \\
            --concentration_parameter ${conc} \\
            ${mcmc_args} \\
            --output    "${gene_id}.FADE.json" \\
        || echo "FADE failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.FADE.json" ] || rm -f "${gene_id}.FADE.json"
        """
    }
}
