#!/usr/bin/env nextflow
// ct_accumulation.nf — Top-level workflow orchestrating CAAS gene accumulation analysis.
// PhyloPhere | workflows/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_ACCUMULATION workflow: aggregates filtered CAAS positions into a gene-level
 *  score matrix, then estimates significance by permutation (randomize phase).
 *
 *  Entry channels (from main.nf):
 *    - params.discovery_dir  : directory of per-trait CT_POSTPROC output TSVs
 *    - params.alignment_dir  : alignment FASTA files (for background position counts)
 *    - params.traitfile      : species × trait matrix
 *
 *  Exit channels (fed to SCORING):
 *    - accumulation_out      : tuple(traitname, direction, path(accumulation_*.tsv))
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

include { CTACC_RUN } from '../subworkflows/CT_ACCUMULATION/ctacc_run.nf'


// ── Workflow ──────────────────────────────────────────────────────────────────

workflow CT_ACCUMULATION {

    take:
    discovery_ch   // tuple(traitname, path(filtered_discovery.tsv))
    alignment_dir  // val(path)
    traitfile      // path

    main:

    // Run accumulation for top, bottom, and all directions in parallel.
    // "all" pools both phenotype extremes, used as a background reference in SCORING.
    directions = Channel.of("top", "bottom", "all")

    accumulation_input = discovery_ch.combine(directions)   // (traitname, disc, direction)

    CTACC_RUN(
        accumulation_input,
        alignment_dir,
        traitfile,
    )

    emit:
    accumulation = CTACC_RUN.out.accumulation   // tuple(traitname, direction, path)
}
