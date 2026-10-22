#!/usr/bin/env nextflow
// rer_cont.nf — RERconverge correlation of gene RERs with a continuous trait.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_CONT: runs continuous_rer.R, which transforms the trait (rer_transform), converts
 *  it to phylogenetic paths on the master tree and correlates the paths with the RER
 *  matrix (correlateWithContinuousPhenotype), with an optional Brownian-motion
 *  permulation null (rer_perm_batches > 0). A failed task stops the workflow.
 *
 *  Consumes:  polished trait RData (trait_vector, n_vector, c_vector), master gene-trees
 *             RDS (RER_TREES), RER matrix RDS (RER_MATRIX)
 *  Produces:  <trait>.char2path.output, <trait>.continuous.output,
 *             <trait>.continuous.perms.rds (directory of corStat blocks, only with permutations)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Continuous correlation ─────────────────────────────────────────────────────

process RER_CONT {
    tag "$rer_matrix"
    errorStrategy 'terminate'
    label 'process_medium'


    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path trait_file
    path rer_master_tree
    path rer_matrix


    output:
    path "${params.traitname}.char2path.output",  emit: char2path
    path "${params.traitname}.continuous.output", emit: continuous_output
    path "${params.traitname}.continuous.perms.rds", emit: perms_output, optional: true



    script:
    def perm_batches    = params.rer_perm_batches      ?: 0
    def perms_per_batch = params.rer_perms_per_batch   ?: 100
    def perm_mode       = 'cc'
    // Include external R sources in the task hash so resume cannot reuse a null
    // computed with a different correlation implementation.
    def sources = ['continuous_rer.R', 'rer_permulations.R', 'rer_batches.R'].collect {
        file("${baseDir}/subworkflows/RERCONVERGE/local/${it}").text
    }.join('\n')
    def source_hash = java.security.MessageDigest.getInstance('SHA-256')
        .digest(sources.getBytes('UTF-8')).encodeHex().toString()

    if (params.use_singularity || params.use_apptainer) {

        """
        echo "Using Singularity/Apptainer"
        # RER sources SHA256: ${source_hash}
        rm -rf -- '${params.traitname}.continuous.perms.rds' '${params.traitname}.continuous.perms.rds.gene_slices'
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/continuous_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.char2path.output \\
        ${rer_matrix} \\
        ${params.traitname}.continuous.output \\
        ${params.rer_minsp} \\
        ${params.winsorizeRER} \\
        ${params.winsorizeTrait} \\
        ${perm_batches} \\
        ${perms_per_batch} \\
        ${perm_mode} \\
        "${params.rer_transform ?: 'auto'}" \\
        ${params.rer_batch_seed}
        """
    } else {
        """
        echo "Running locally"
        # RER sources SHA256: ${source_hash}
        rm -rf -- '${params.traitname}.continuous.perms.rds' '${params.traitname}.continuous.perms.rds.gene_slices'
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/continuous_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.char2path.output \\
        ${rer_matrix} \\
        ${params.traitname}.continuous.output \\
        ${params.rer_minsp} \\
        ${params.winsorizeRER} \\
        ${params.winsorizeTrait} \\
        ${perm_batches} \\
        ${perms_per_batch} \\
        ${perm_mode} \\
        "${params.rer_transform ?: 'auto'}" \\
        ${params.rer_batch_seed}
        """
    }

}
