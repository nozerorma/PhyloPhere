#!/usr/bin/env nextflow
// rer_matrix.nf — Matrix of relative evolutionary rates (RERs) of all genes.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_MATRIX: runs rer_matrix.R, which computes the RER of every gene on every branch
 *  of the master tree (getAllResiduals) for the species of the trait vector.
 *
 *  Consumes:  polished trait RData (trait_vector), master gene-trees RDS (RER_TREES)
 *  Produces:  <trait>.RERmatrix.output (RDS)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── RER matrix ─────────────────────────────────────────────────────────────────

process RER_MATRIX {
    tag "$gene_trees_file"

    label 'process_medium'


    publishDir path: "${params.outdir}/rerconverge/rer_objects", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path trait_file
    path gene_trees_file

    output:
    file("${params.traitname}.RERmatrix.output")


    script:
    // Extra arguments from task.ext.args
    def args = task.ext.args ?: ''
    def matrix_out = "${params.traitname}.RERmatrix.output"

    if (params.use_singularity) {
        """
        echo "Using Singularity"
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/rer_matrix.R' \\
        ${ trait_file } \\
        ${ gene_trees_file } \\
        ${ matrix_out } \\
        ${params.rer_minsp} \\
        $args
        """
    } else {
        """
        echo "Running locally"
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/rer_matrix.R' \\
        ${ trait_file } \\
        ${ gene_trees_file } \\
        ${ matrix_out } \\
        ${params.rer_minsp} \\
        $args
        """
    }

}
