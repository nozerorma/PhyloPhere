#!/usr/bin/env nextflow
// rer_bin.nf — RERconverge correlation of gene RERs with a binary (0/1) trait.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_BIN: runs binary_rer.R, which builds the foreground paths of the 0/1 trait on
 *  the master tree and correlates them with the RER matrix
 *  (correlateWithBinaryPhenotype), with an optional permulation null
 *  (rer_perm_batches > 0). A failed task stops the workflow.
 *
 *  Consumes:  polished trait RData (trait_vector, 0/1), master gene-trees RDS (RER_TREES),
 *             RER matrix RDS (RER_MATRIX)
 *  Produces:  <trait>.fg_paths.output, <trait>.binary.output,
 *             <trait>.binary.perms.rds (directory of corStat blocks, only with permutations)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Binary correlation ─────────────────────────────────────────────────────────

process RER_BIN {
    tag "$rer_matrix"
    label 'process_medium'
    errorStrategy 'terminate'

    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path trait_file
    path rer_master_tree
    path rer_matrix

    output:
    path "${params.traitname}.fg_paths.output",  emit: fg_paths
    path "${params.traitname}.binary.output",    emit: binary_output
    path "${params.traitname}.binary.perms.rds",  emit: perms_output, optional: true

    script:
    def perm_batches    = params.rer_perm_batches    ?: 0
    def perms_per_batch = params.rer_perms_per_batch ?: 100
    def min_pos         = params.rer_min_pos == null ? 2 : params.rer_min_pos
    def binary_clade    = params.rer_binary_clade    ?: 'all'
    // Include external R sources in the task hash so resume cannot reuse a null
    // computed with a different correlation implementation.
    def sources = ['binary_rer.R', 'rer_permulations.R', 'rer_batches.R'].collect {
        file("${baseDir}/subworkflows/RERCONVERGE/local/${it}").text
    }.join('\n')
    def source_hash = java.security.MessageDigest.getInstance('SHA-256')
        .digest(sources.getBytes('UTF-8')).encodeHex().toString()

    if (params.use_singularity || params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        # RER sources SHA256: ${source_hash}
        rm -rf -- '${params.traitname}.binary.perms.rds' '${params.traitname}.binary.perms.rds.gene_slices'
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/binary_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.fg_paths.output \\
        ${rer_matrix} \\
        ${params.traitname}.binary.output \\
        ${params.rer_minsp} \\
        ${min_pos} \\
        ${params.winsorizeRER} \\
        ${binary_clade} \\
        ${perm_batches} \\
        ${perms_per_batch} \\
        ${params.rer_batch_seed}
        """
    } else {
        """
        echo "Running locally"
        # RER sources SHA256: ${source_hash}
        rm -rf -- '${params.traitname}.binary.perms.rds' '${params.traitname}.binary.perms.rds.gene_slices'
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/binary_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.fg_paths.output \\
        ${rer_matrix} \\
        ${params.traitname}.binary.output \\
        ${params.rer_minsp} \\
        ${min_pos} \\
        ${params.winsorizeRER} \\
        ${binary_clade} \\
        ${perm_batches} \\
        ${perms_per_batch} \\
        ${params.rer_batch_seed}
        """
    }
}
