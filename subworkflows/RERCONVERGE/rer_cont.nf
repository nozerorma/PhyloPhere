#!/usr/bin/env nextflow
// rer_cont.nf — RERconverge correlation of gene RERs with a continuous trait.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_CONT: runs continuous_rer.R, which transforms the trait (rer_transform), converts
 *  it to phylogenetic paths on the master tree and correlates the paths with the RER
 *  matrix (correlateWithContinuousPhenotype), with an optional Brownian-motion
 *  permulation null (rer_perm_batches > 0). A failed task is ignored.
 *
 *  Consumes:  polished trait RData (trait_vector, n_vector, c_vector), master gene-trees
 *             RDS (RER_TREES), RER matrix RDS (RER_MATRIX)
 *  Produces:  <trait>.char2path.output, <trait>.continuous.output,
 *             <trait>.continuous.perms.rds (only with permutations)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Continuous correlation ─────────────────────────────────────────────────────

process RER_CONT {
    tag "$rer_matrix"
    errorStrategy 'ignore'
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

    if (params.use_singularity || params.use_apptainer) {

        """
        echo "Using Singularity/Apptainer"
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
        "${params.rer_transform ?: 'auto'}"
        """
    } else {
        """
        echo "Running locally"
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
        "${params.rer_transform ?: 'auto'}"
        """
    }

}
