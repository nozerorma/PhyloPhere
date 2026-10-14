#!/usr/bin/env nextflow
// fade_tables.nf — The published FADE tables of one direction, merged from the tables of the FADE tasks.
// PhyloPhere | subworkflows/FADE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  FADE_TABLES: concatenates the gene and site tables of every FADE task (FADE_BATCHED, FADE_RUN) and keeps the rows of one
 *  direction (top or bottom), through fade_merge_tables.py.
 *
 *  Consumes:  direction, the collected *.fade_genes.tsv and *.fade_sites.tsv files (empty when no task ran)
 *  Produces:  fade_summary_<direction>.tsv      one row per gene (read by SCORING and FADE_GENE_LISTS)
 *             fade_gene_counts_<direction>.tsv  the same rows with the sigcount_<AA> columns (read by FADE_REPORT)
 *             fade_site_bf_<direction>.tsv      one row per significant site: gene, position (0-based), max_bf, target_aa (read by SCORING)
 *             fade_sites_<direction>.csv        the same rows, comma-delimited (read by SCORING and ENRICHMENT)
 *  Every table has its header when no gene reached it.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

process FADE_TABLES {
    tag "fade_tables|${direction}"
    label 'process_low'

    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_summary_*.tsv'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_gene_counts_*.tsv'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_site_bf_*.tsv'
    publishDir path: "${params.outdir}/selection/fade/${direction}",
               mode: 'copy', overwrite: true,
               pattern: 'fade_sites_*.csv'

    input:
    val  direction
    path genes
    path sites

    output:
    path "fade_summary_${direction}.tsv",      emit: summary_tsv
    path "fade_gene_counts_${direction}.tsv",  emit: gene_counts_tsv
    path "fade_site_bf_${direction}.tsv",      emit: site_tsv
    path "fade_sites_${direction}.csv",        emit: sites_csv

    script:
    def pyBin = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    """
    ${pyBin} ${baseDir}/subworkflows/FADE/local/src/fade_merge_tables.py \\
        --direction ${direction} --genes ${genes} --sites ${sites}
    """
}
