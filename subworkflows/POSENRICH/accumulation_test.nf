#!/usr/bin/env nextflow
// accumulation_test.nf — Exact test of CAAS accumulation per gene and direction.
// PhyloPhere | subworkflows/POSENRICH/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_TEST: runs local/src/accumulation_test.py once (level 1, share of detections per gene
 *  calibrated with the null cycles). The script block has a container branch (entrypoint wrapper) and a plain
 *  branch that run the same command.
 *
 *  Consumes:  universe.tsv.gz of POSENRICH_UNIVERSE (the universe U shared with level 2) and the block
 *             directories of the null (positions.parquet, cycles.txt)
 *  Produces:  accumulation/accumulation_results.tsv, accumulation/accumulation_manifest.json
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


process ACCUMULATION_TEST {
    tag "accumulation_test"
    label 'process_medium'

    publishDir path: { "${params.outdir}/posenrich" }, mode: 'copy', overwrite: true, pattern: 'accumulation'

    input:
    path universe               // universe.tsv.gz: one row per position of U, with the observed positive flags
    path blocks, stageAs: 'blk_*'   // block directories of the null (positions.parquet, cycles.txt)

    output:
    path 'accumulation',        emit: results

    script:
    def script  = "${baseDir}/subworkflows/POSENRICH/local/src/accumulation_test.py"
    def wrapper = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh ' : ''
    """
    ${wrapper}python ${script} \\
        --universe '${universe}' --blocks blk_* --output-dir accumulation
    """
}
