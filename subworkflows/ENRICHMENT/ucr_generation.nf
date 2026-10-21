#!/usr/bin/env nextflow
// ucr_generation.nf — Generate ucr_positions.tsv (ultra-conserved regions) from the protein alignments.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  UCR_GENERATION: builds the ucr_positions.tsv that POSENRICH takes as
 *  --ucr_positions_file when the file is not provided. Three steps:
 *    1. COMPUTE_ALIGNMENT_ENTROPY (compute_alignment_entropy.py): per-gene Valdar
 *       variability.
 *    2. RUN_UCR_DETECTION (run_ucr_detection.py, which calls detect_ucr.py): UCR
 *       windows per gene.
 *    3. AGGREGATE_UCR_POSITIONS (aggregate_ucr.py): ucr_positions.tsv with the columns
 *       build_position_gmt.py reads (gene, ucr_id, method, position, region_type,
 *       C_trident, variability, g).
 *
 *  The caller runs it only when --tax_id and --alignment are available (the entropy
 *  step needs the taxonomy table); otherwise no UCR layer is built, which POSENRICH
 *  accepts as an absent optional layer.
 *
 *  Consumes:  alignment directory, tax_id table
 *  Produces:  core_inputs/ucr_positions.tsv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Alignment variability ──────────────────────────────────────────────────────

// Valdar variability per column and gene (compute_alignment_entropy.py)
process COMPUTE_ALIGNMENT_ENTROPY {
    tag "auto-generate Valdar variability"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path alignment_dir
    path taxid_tsv

    output:
    path "entropy_dir", emit: entropy_dir

    script:
    """
    python3 ${baseDir}/bin/compute_alignment_entropy.py \\
        --alignment-dir "${alignment_dir}" \\
        --output-dir entropy_dir \\
        --taxid-tsv "${taxid_tsv}"
    """
}

// ── UCR detection ──────────────────────────────────────────────────────────────

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


// ── Aggregation ────────────────────────────────────────────────────────────────

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


// ── Workflow ───────────────────────────────────────────────────────────────────

workflow UCR_GENERATION {
    take:
        alignment_dir_ch  // alignment directory
        taxid_ch          // tax_id table

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
