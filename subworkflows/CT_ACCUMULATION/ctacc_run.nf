#!/usr/bin/env nextflow
// ctacc_run.nf — Exact test of CAAS accumulation per gene and direction.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_ACCUMULATION_TEST: runs local/accumulation_test.py once per direction (top, bottom, all).
 *  The script block has a container branch (entrypoint wrapper) and a plain branch that run
 *  the same command.
 *
 *  Consumes:  filtered_discovery.tsv and cleaned background list (CT_POSTPROC), the tested
 *             positions (background.output) and null_positions.parquet of one
 *             CAAS_CORE_MERGE_REDUCE run (with the observed positions, the universe U)
 *  Produces:  accumulation_<direction>_results.csv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// One run per direction ('top', 'bottom', 'all'): the exact hypergeometric test of every gene on the universe U
process CT_ACCUMULATION_TEST {
    tag "ct_accumulation_test|${direction}"
    label 'process_medium'

    publishDir path: { "${params.outdir}/accumulation/${direction}" },
               mode: 'copy', overwrite: true,
               pattern: '*_results.csv'

    input:
    val  direction
    path caas_csv
    path background_positions   // caastools background.output (tested positions)
    path bg_caas_universe       // cleaned_background_main.txt (surviving genes)
    path null_positions         // null_positions.parquet from CAAS_CORE_MERGE_REDUCE, the filtered null SCORING also uses:
                                 // the tested positions that score in it, plus the observed ones, are the universe U

    output:
    val  direction,                           emit: direction
    path "accumulation_*_results.csv",        emit: results

    script:
    def script       = "${baseDir}/subworkflows/CT_ACCUMULATION/local/accumulation_test.py"
    def out_file     = "accumulation_${direction}_results.csv"
    // 'all' maps to --change-side both: a position detected on both sides counts once
    def change_side_arg = (direction == 'all') ? 'both' : direction
    def wrapper      = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh ' : ''
    """
    ${wrapper}python ${script} \\
        --caas-csv '${caas_csv}' \\
        --background-positions '${background_positions}' \\
        --bg-caas '${bg_caas_universe}' \\
        --null-positions '${null_positions}' \\
        --change-side '${change_side_arg}' \\
        --output '${out_file}'
    """
}
