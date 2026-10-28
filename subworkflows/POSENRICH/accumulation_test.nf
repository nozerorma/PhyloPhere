#!/usr/bin/env nextflow
// accumulation_test.nf — Exact test of CAAS accumulation per gene and direction.
// PhyloPhere | subworkflows/POSENRICH/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_TEST: runs local/src/accumulation_test.py once per direction (top, bottom, all).
 *  The script block has a container branch (entrypoint wrapper) and a plain branch that run
 *  the same command.
 *
 *  Consumes:  universe.tsv.gz of POSENRICH_UNIVERSE (the universe U shared with level 2)
 *  Produces:  accumulation_<direction>_results.csv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// One run per direction ('top', 'bottom', 'all'): the exact hypergeometric test of every gene of the universe U
process ACCUMULATION_TEST {
    tag "accumulation_test|${direction}"
    label 'process_medium'

    publishDir path: { "${params.outdir}/posenrich/accumulation/${direction}" },
               mode: 'copy', overwrite: true,
               pattern: '*_results.csv'

    input:
    val  direction
    path universe               // universe.tsv.gz: one row per position of U, with the observed positive flags

    output:
    val  direction,                           emit: direction
    path "accumulation_*_results.csv",        emit: results

    script:
    def script       = "${baseDir}/subworkflows/POSENRICH/local/src/accumulation_test.py"
    def out_file     = "accumulation_${direction}_results.csv"
    // 'all' maps to --change-side both: a position detected on both sides counts once
    def change_side_arg = (direction == 'all') ? 'both' : direction
    def wrapper      = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh ' : ''
    """
    ${wrapper}python ${script} \\
        --universe '${universe}' \\
        --change-side '${change_side_arg}' \\
        --output '${out_file}'
    """
}
