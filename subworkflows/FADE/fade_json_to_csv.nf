#!/usr/bin/env nextflow
// fade_json_to_csv.nf — Site-level FADE table (gene, position, max_bf, target_aa) from the raw JSON files.
// PhyloPhere | subworkflows/FADE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  FADE_JSON_TO_CSV: parses the *.FADE.json files of one direction into one CSV with a
 *  row per site whose maximum Bayes factor over the target amino acids reaches
 *  fade_bf_threshold. The gene-level FADE report keeps no site-level table, so this
 *  is the position-keyed ((Gene, Position)) FADE evidence that POSENRICH uses, like
 *  the UCR and FUBAR layers. It is not gated behind --enrichment.
 *
 *  Consumes:  direction ('top' or 'bottom'), the collected *.FADE.json files
 *  Produces:  fade_sites_<direction>.csv (gene, position, max_bf, target_aa; the
 *             header row is written even when no site reaches the threshold)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── JSON to site table ─────────────────────────────────────────────────────────

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
