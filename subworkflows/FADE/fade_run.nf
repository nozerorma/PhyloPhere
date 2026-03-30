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

process FADE_BATCHED {
    tag "$batchID (${batchSize} genes, ${direction})"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/fade/${direction}/json",
               mode: 'copy', overwrite: true,
               pattern: '*.FADE.json'

    // Batch-level success is independent of individual gene failures —
    // the batch script logs failures and continues.

    input:
    tuple val(batchID), val(direction), val(batchSize), val(batchManifestText),
          path(fastas, stageAs: 'fastas/*'), path(trees, stageAs: 'trees/*')
    path lg_dat

    output:
    tuple val(direction), path("*.FADE.json"), emit: fade_json, optional: true

    script:
    def model        = params.fade_model        ?: 'LG'
    def model_file_arg = model == 'LG' ? "--model-file lg.dat" : ""
    def method       = params.fade_method       ?: 'Variational-Bayes'
    def grid         = params.fade_grid         ?: 20
    def conc         = params.fade_concentration ?: 0.5
    def runnerMode   = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'

    def mcmc_args = (method == 'Variational-Bayes') ? "" :
        """--mcmc-chains ${params.fade_chains ?: 5} \\
           --mcmc-chain-length ${params.fade_chain_length ?: 2000000} \\
           --mcmc-burn-in ${params.fade_burn_in ?: 1000000} \\
           --mcmc-samples ${params.fade_samples ?: 1000}"""

    """
cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

bash ${baseDir}/subworkflows/FADE/local/scripts/run_hyphy_fade_batch.sh \\
    --batch-id       ${batchID} \\
    --manifest       ${batchID}.manifest.tsv \\
    --direction      ${direction} \\
    --workers        ${params.fade_batch_workers} \\
    --runner-mode    ${runnerMode} \\
    --model          ${model} \\
    --model-file-arg "${model_file_arg}" \\
    --method         "${method}" \\
    --grid           ${grid} \\
    --concentration  ${conc} \\
    ${mcmc_args}
"""
}

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
    tuple val(gene_id), val(direction), path("${gene_id}.${direction}.FADE.json"), emit: fade_json, optional: true

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
            --output    "${gene_id}.${direction}.FADE.json" \\
        || echo "FADE failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.FADE.json" ] || rm -f "${gene_id}.${direction}.FADE.json"
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
            --output    "${gene_id}.${direction}.FADE.json" \\
        || echo "FADE failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.FADE.json" ] || rm -f "${gene_id}.${direction}.FADE.json"
        """
    }
}
