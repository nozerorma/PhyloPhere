#!/usr/bin/env nextflow

/*
 * MOLERATE_RUN
 * ────────────
 * Run HyPhy MoleRate (relative evolutionary molecular rate) on a single gene
 * alignment, comparing evolutionary rates between foreground and background
 * branches relative to a shared reference tree with branch lengths.
 *
 * MoleRate fits four nested models (Proportional, Proportional Partitioned,
 * Unconstrained Test, Unconstrained) and reports likelihood-ratio tests
 * comparing them, along with per-test-branch rate estimates.
 *
 * Inputs
 * ──────
 *   gene_id    : string  — gene identifier
 *   direction  : string  — 'top' or 'bottom'
 *   fasta      : protein FASTA alignment
 *   tree       : reference Newick tree with branch lengths
 *   fg_list    : plain-text file with FG species names (one per line)
 *
 * Outputs
 * ───────
 *   molerate_json : MoleRate results JSON
 */

process MOLERATE_BATCHED {
    tag "$batchID (${batchSize} genes, ${direction})"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/molerate/${direction}/json",
               mode: 'copy', overwrite: true,
               pattern: '*.molerate.json'

    // Batch-level success is independent of individual gene failures —
    // the batch script logs failures and continues.

    input:
    tuple val(batchID), val(direction), val(batchSize), val(batchManifestText),
          path(fastas, stageAs: 'fastas/*'), path(tree), path(fg_list)
    path lg_dat

    output:
    tuple val(direction), path("*.molerate.json"), emit: molerate_json, optional: true

    script:
    def hyphy_analyses = params.hyphy_analyses_path ?: '/hyphy-analyses'
    def model          = params.molerate_model        ?: 'LG'
    def model_file_arg = model == 'LG' ? "--model-file lg.dat" : ""
    def rv             = params.molerate_rv           ?: 'None'
    def rate_classes   = params.molerate_rate_classes ?: 4
    def labeling       = params.molerate_labeling_strategy ?: 'all-descendants'
    def branch_tests   = (params.molerate_branch_tests != false) ? 'Yes' : 'No'
    def full_model     = (params.molerate_skip_full_model == true) ? 'No' : 'Yes'
    def runnerMode     = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'

    """
cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

bash ${baseDir}/subworkflows/MOLERATE/local/scripts/run_hyphy_molerate_batch.sh \\
    --batch-id           ${batchID} \\
    --manifest           ${batchID}.manifest.tsv \\
    --direction          ${direction} \\
    --workers            ${params.molerate_batch_workers} \\
    --runner-mode        ${runnerMode} \\
    --tree               "${tree}" \\
    --fg-list            "${fg_list}" \\
    --hyphy-analyses     ${hyphy_analyses} \\
    --model              ${model} \\
    --model-file-arg     "${model_file_arg}" \\
    --rv                 ${rv} \\
    --rate-classes       ${rate_classes} \\
    --labeling-strategy  ${labeling} \\
    --branch-tests       ${branch_tests} \\
    --full-model         ${full_model}
"""
}

process MOLERATE_RUN {
    tag "${gene_id}|${direction}"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/molerate/${direction}/json",
               mode: 'copy', overwrite: true,
               pattern: '*.molerate.json'

    errorStrategy 'ignore'  // Skip genes that fail (e.g., insufficient branch resolution)

    input:
    tuple val(gene_id), val(direction), path(fasta), path(tree), path(fg_list)
    path lg_dat

    output:
    tuple val(gene_id), val(direction), path("${gene_id}.${direction}.molerate.json"), emit: molerate_json, optional: true

    script:
    def hyphy_analyses = params.hyphy_analyses_path ?: '/hyphy-analyses'
    def model          = params.molerate_model        ?: 'LG'
    def model_file_arg = model == 'LG' ? "--model-file lg.dat" : ""
    def rv             = params.molerate_rv           ?: 'None'
    def rate_classes   = params.molerate_rate_classes ?: 4
    def labeling       = params.molerate_labeling_strategy ?: 'all-descendants'
    def branch_tests   = (params.molerate_branch_tests != false) ? 'Yes' : 'No'
    def full_model     = (params.molerate_skip_full_model == true) ? 'No' : 'Yes'

    // Build repeated --branches args from the fg_list file at runtime
    if (params.use_singularity || params.use_apptainer) {
        """
        # Build --branches arguments from fg_list
        BRANCHES_ARGS=\$(awk '{printf " --branches %s", \$1}' "${fg_list}")

        /usr/local/bin/_entrypoint.sh hyphy "${hyphy_analyses}/molerate/molerate.bf" \\
            --alignment "${fasta}" \\
            --tree      "${tree}" \\
            --model     ${model} \\
            ${model_file_arg} \\
            --rv        ${rv} \\
            --rate-classes ${rate_classes} \\
            --labeling-strategy ${labeling} \\
            --branch-level-analysis ${branch_tests} \\
            --full-model ${full_model} \\
            --output "${gene_id}.${direction}.molerate.json" \\
            \$BRANCHES_ARGS \\
        || echo "MoleRate failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.molerate.json" ] || rm -f "${gene_id}.${direction}.molerate.json"
        """
    } else {
        """
        BRANCHES_ARGS=\$(awk '{printf " --branches %s", \$1}' "${fg_list}")

        # Resolve molerate.bf: prefer configured path; fall back to HyPhy's own
        # TemplateBatchFiles when running locally without a container.
        MOLERATE_BF="${hyphy_analyses}/molerate/molerate.bf"
        if [ ! -f "\$MOLERATE_BF" ]; then
            MOLERATE_BF="\$(dirname \$(which hyphy))/../share/hyphy/TemplateBatchFiles/molerate.bf"
        fi

        hyphy "\$MOLERATE_BF" \\
            --alignment "${fasta}" \\
            --tree      "${tree}" \\
            --model     ${model} \\
            ${model_file_arg} \\
            --rv        ${rv} \\
            --rate-classes ${rate_classes} \\
            --labeling-strategy ${labeling} \\
            --branch-level-analysis ${branch_tests} \\
            --full-model ${full_model} \\
            --output "${gene_id}.${direction}.molerate.json" \\
            \$BRANCHES_ARGS \\
        || echo "MoleRate failed for ${gene_id} (${direction}), skipping"
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.molerate.json" ] || rm -f "${gene_id}.${direction}.molerate.json"
        """
    }
}
