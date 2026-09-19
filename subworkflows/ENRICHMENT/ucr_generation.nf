#!/usr/bin/env nextflow

/*
 * UCR_GENERATION subworkflow
 *
 * Auto-generates ucr_positions.tsv (posenrich's --ucr_positions_file) when
 * left blank, chaining three verbatim ports of the user's own reference
 * implementation (ortholog_characterizator/subworkflows/variability):
 *   1. compute_alignment_entropy.py -> per-gene Valdar variability
 *      (shared with CT_ACCUMULATION's own auto-generation, §5)
 *   2. run_ucr_detection.py -> bin/detect_ucr.py (verbatim) per gene
 *   3. bin/aggregate_ucr.py (verbatim) -> ucr_positions.tsv in the exact
 *      schema build_position_gmt.py expects (gene, ucr_id, method, position,
 *      region_type, C_trident, variability, g)
 *
 * Requires --tax_id (compute_variability.py's hard requirement) and
 * --alignment. If either is unavailable, generation is skipped and
 * ucr_positions_file stays empty — posenrich_enrich.py already treats an
 * absent UCR file as an optional annotation layer, not a hard failure.
 *
 * NOTE: if --ct_accumulation is ALSO enabled in the same run, its own
 * COMPUTE_ALIGNMENT_ENTROPY call (workflows/ct_accumulation.nf) recomputes
 * the identical entropy files independently — not deduplicated across
 * workflows in a fresh run (Nextflow's own -resume work-dir cache does
 * dedupe it on a resumed run, since the process+inputs hash is identical).
 * Not worth a cross-module cache for a single redundant computation.
 */

include { COMPUTE_ALIGNMENT_ENTROPY } from "${baseDir}/subworkflows/CT_ACCUMULATION/ctacc_run.nf"

process RUN_UCR_DETECTION {
    tag "auto-generate UCR windows"
    label 'process_medium'

    input:
    path entropy_dir

    output:
    path "ucr_raw", emit: ucr_raw_dir

    script:
    """
    python3 ${baseDir}/bin/run_ucr_detection.py \\
        --entropy-dir "${entropy_dir}" \\
        --output-dir ucr_raw
    """
}

process AGGREGATE_UCR_POSITIONS {
    tag "auto-generate ucr_positions.tsv"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path ucr_raw_dir
    path entropy_dir
    path alignment_dir
    path taxid_tsv

    output:
    path "ucr_agg/ucr_positions.tsv", emit: ucr_positions_file

    script:
    """
    mkdir -p ucr_agg
    python3 ${baseDir}/bin/aggregate_ucr.py \\
        --ucr_dir "${ucr_raw_dir}" \\
        --entropy_dir "${entropy_dir}" \\
        --prot_dir "${alignment_dir}" \\
        --taxid_tsv "${taxid_tsv}" \\
        --out_dir ucr_agg
    """
}

workflow UCR_GENERATION {
    take:
        alignment_dir_ch  // path/value channel: alignment directory
        taxid_ch          // path/value channel: tax_id file

    main:
        COMPUTE_ALIGNMENT_ENTROPY(alignment_dir_ch, taxid_ch)
        RUN_UCR_DETECTION(COMPUTE_ALIGNMENT_ENTROPY.out.entropy_dir)
        AGGREGATE_UCR_POSITIONS(
            RUN_UCR_DETECTION.out.ucr_raw_dir,
            COMPUTE_ALIGNMENT_ENTROPY.out.entropy_dir,
            alignment_dir_ch,
            taxid_ch,
        )

    emit:
        ucr_positions_file = AGGREGATE_UCR_POSITIONS.out.ucr_positions_file
}
