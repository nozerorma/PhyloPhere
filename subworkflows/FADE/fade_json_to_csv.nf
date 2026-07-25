#!/usr/bin/env nextflow

/*
 * FADE_JSON_TO_CSV
 * ─────────────────
 * Parses the raw *.FADE.json output directly into a simple per-site CSV
 * (gene, position, max_bf, target_aa) for one direction. HyPhy FADE's
 * gene-level report no longer keeps site-level detail (dropped as a memory
 * optimisation once nothing downstream consumed it) — this restores exactly
 * the piece posenrich needs: a position-keyed ((Gene,Position)) evidence
 * layer, analogous to UCR core/flank and FUBAR positive/purifying, built
 * from the classic FADE significance criterion (BF >= fade_bf_threshold).
 *
 * Runs unconditionally (not gated behind --enrichment/--string) so the
 * position-level FADE evidence is always available for posenrich.
 *
 * Inputs
 * ──────
 *   direction   : val  — 'top' or 'bottom'
 *   json_files  : collected list of *.FADE.json paths
 *
 * Outputs
 * ───────
 *   sites_csv   : gene,position,max_bf,target_aa (header row always present,
 *                 even when zero sites clear the threshold)
 */

process FADE_JSON_TO_CSV {
    tag "fade_json_to_csv|${direction}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir path: { "${params.outdir}/selection/fade/${direction}" },
               mode: 'copy', overwrite: true,
               pattern: 'fade_sites_*.csv'

    input:
    val  direction
    path json_files

    output:
    path "fade_sites_${direction}.csv", emit: sites_csv

    script:
    def bf_thr  = params.fade_bf_threshold ?: 100
    def n_cores = task.cpus ?: 4
    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh Rscript ${baseDir}/subworkflows/FADE/local/src/parse_fade_json_sites.R \
            --json_dir  . \
            --direction ${direction} \
            --bf_thr    ${bf_thr} \
            --n_cores   ${n_cores} \
            --out       fade_sites_${direction}.csv
        """
    } else {
        """
        Rscript ${baseDir}/subworkflows/FADE/local/src/parse_fade_json_sites.R \
            --json_dir  . \
            --direction ${direction} \
            --bf_thr    ${bf_thr} \
            --n_cores   ${n_cores} \
            --out       fade_sites_${direction}.csv
        """
    }
}
