#!/usr/bin/env nextflow
// fade_run.nf — HyPhy FADE (directional amino-acid selection) on foreground-annotated gene trees.
// PhyloPhere | subworkflows/FADE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  FADE_RUN, FADE_BATCHED: run HyPhy FADE (FUBAR Approach to Directional Evolution) on
 *  a protein alignment and a tree whose foreground branches are labeled {Foreground}
 *  (from ANNOTATE_TREE_FG), testing the foreground branches for directional selection
 *  toward each amino acid. FADE_RUN handles one gene; FADE_BATCHED handles
 *  fade_batch_size genes per task through run_hyphy_fade_batch.sh. A gene whose FADE
 *  run fails is skipped and produces no JSON.
 *
 *  Consumes:  gene id, direction ('top' or 'bottom'), protein alignment (taxa names
 *             equal to the tree labels), annotated tree, lg_dat (model .dat file
 *             staged in the task directory)
 *  Produces:  <gene>.<direction>.FADE.json (HyPhy standard output)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Batched run ────────────────────────────────────────────────────────────────

process FADE_BATCHED {
    tag "$batchID (${batchSize} genes, ${direction})"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/fade/top/json",
               mode: 'copy', overwrite: true,
               pattern: '*.top.FADE.json'
    publishDir path: "${params.outdir}/selection/fade/bottom/json",
               mode: 'copy', overwrite: true,
               pattern: '*.bottom.FADE.json'

    // A failed gene does not fail the task: the batch script logs it and continues.

    input:
    tuple val(batchID), val(direction), val(batchSize), val(batchManifestText),
          path(fastas, stageAs: 'fastas/*'), path(trees, stageAs: 'trees/*')
    path lg_dat

    output:
    path "*.FADE.json", emit: fade_json, optional: true

    script:
    def model         = params.fade_model        ?: 'GTR'
    def method        = params.fade_method       ?: 'Variational-Bayes'
    def grid          = params.fade_grid         ?: 20
    def conc          = params.fade_concentration ?: 0.5
    def shareBaseline = ((params.fade_background_scope ?: 'all') == 'all') ? 'true' : 'false'
    def runnerMode    = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'
    def nWorkers      = (task.cpus ?: 8) as int
    // Each HyPhy call gets floor(task.cpus / workers) CPUs, at least 1, so that it does
    // not read the CPU count of the whole node and oversubscribe it.
    def cpuPerWorker  = Math.max(1, (task.cpus as int).intdiv(nWorkers))

    def mcmc_args = (method == 'Variational-Bayes') ? "" :
        """--mcmc-chains ${params.fade_chains ?: 5} \\
           --mcmc-chain-length ${params.fade_chain_length ?: 2000000} \\
           --mcmc-burn-in ${params.fade_burn_in ?: 1000000} \\
           --mcmc-samples ${params.fade_samples ?: 100}"""

    """
cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

bash ${baseDir}/subworkflows/FADE/local/src/run_hyphy_fade_batch.sh \\
    --batch-id       ${batchID} \\
    --manifest       ${batchID}.manifest.tsv \\
    --direction      ${direction} \\
    --workers        ${nWorkers} \\
    --cpu-per-worker ${cpuPerWorker} \\
    --runner-mode    ${runnerMode} \\
    --model          ${model} \\
    --method         "${method}" \\
    --grid           ${grid} \\
    --concentration  ${conc} \\
    --share-baseline ${shareBaseline} \\
    ${mcmc_args}
"""
}


// ── Single-gene run ────────────────────────────────────────────────────────────

process FADE_RUN {
    tag "${gene_id}|${direction}"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/fade/${direction}/json",
               mode: 'copy', overwrite: true,
               pattern: '*.FADE.json'

    errorStrategy 'ignore'  // a gene that fails (e.g. too few foreground branches) is skipped

    input:
    tuple val(gene_id), val(direction), path(fasta), path(annotated_tree)
    path lg_dat

    output:
    tuple val(gene_id), val(direction), path("${gene_id}.${direction}.FADE.json"), emit: fade_json, optional: true

    script:
    def model   = params.fade_model   ?: 'GTR'
    def method  = params.fade_method  ?: 'Variational-Bayes'
    def grid    = params.fade_grid    ?: 20
    def conc    = params.fade_concentration ?: 0.5
    def fade_bf = "${baseDir}/subworkflows/FADE/local/src/FADE_phylophere.bf"

    // MCMC options apply only when the method is not Variational-Bayes
    def mcmc_args = (method == 'Variational-Bayes') ? "" :
        """--chains ${params.fade_chains ?: 5} \\
           --chain-length ${params.fade_chain_length ?: 2000000} \\
           --burn-in ${params.fade_burn_in ?: 1000000} \\
           --samples ${params.fade_samples ?: 100}"""

    if (params.use_singularity || params.use_apptainer) {
        """
        export OMP_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}
        export BLAS_NUM_THREADS=${task.cpus}
        /usr/local/bin/_entrypoint.sh hyphy "${fade_bf}" \\
            --alignment "${fasta}" \\
            --tree      "${annotated_tree}" \\
            --branches  Foreground \\
            --model     ${model} \\
            --method    "${method}" \\
            --grid      ${grid} \\
            --concentration_parameter ${conc} \\
            --cache     NONE \\
            --cpu       ${task.cpus} \\
            ${mcmc_args} \\
            --output    "${gene_id}.${direction}.FADE.json" > "${gene_id}.${direction}.hyphy.log" 2>&1 \\
        || { tail -n 30 "${gene_id}.${direction}.hyphy.log" >&2 || true; echo "FADE failed for ${gene_id} (${direction}), skipping"; }
        rm -f "${gene_id}.${direction}.hyphy.log" *.cache
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.FADE.json" ] || rm -f "${gene_id}.${direction}.FADE.json"
        """
    } else {
        """
        export OMP_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}
        export BLAS_NUM_THREADS=${task.cpus}
        hyphy "${fade_bf}" \\
            --alignment "${fasta}" \\
            --tree      "${annotated_tree}" \\
            --branches  Foreground \\
            --model     ${model} \\
            --method    "${method}" \\
            --grid      ${grid} \\
            --concentration_parameter ${conc} \\
            --cache     NONE \\
            --cpu       ${task.cpus} \\
            ${mcmc_args} \\
            --output    "${gene_id}.${direction}.FADE.json" > "${gene_id}.${direction}.hyphy.log" 2>&1 \\
        || { tail -n 30 "${gene_id}.${direction}.hyphy.log" >&2 || true; echo "FADE failed for ${gene_id} (${direction}), skipping"; }
        rm -f "${gene_id}.${direction}.hyphy.log" *.cache
        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "${gene_id}.${direction}.FADE.json" ] || rm -f "${gene_id}.${direction}.FADE.json"
        """
    }
}
