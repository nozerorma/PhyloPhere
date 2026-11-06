#!/usr/bin/env nextflow
include { selectionQuote; selectionHash } from '../SELECTION/selection_contract'
include { fadeJson } from './fade_contract'

// Space-separated, quoted file list for the command line (a function, not a closure, for the linter).
def fadeArgs(paths) {
    return (paths instanceof List ? paths : [paths]).collect { selectionQuote(it) }.join(' ')
}
process FADE_TABLES {
    tag "${direction}"
    label 'process_low'
    cache true
    publishDir "${params.outdir}/selection/fade/${direction}", mode:'copy', overwrite:true
    input:
    val direction
    path genes, stageAs:'genes/*'
    path sites, stageAs:'sites/*'
    path audits, stageAs:'audits/*'
    val expected
    val settings
    val source_hash
    output:
    path "fade_summary_${direction}.tsv", emit: summary_tsv
    path "fade_gene_counts_${direction}.tsv", emit: gene_counts_tsv
    path "fade_site_bf_${direction}.tsv", emit: site_tsv
    path "fade_sites_${direction}.csv", emit: sites_csv
    path 'fade_run_metadata.json', emit: metadata
    path 'fade_audit.tsv', emit: audit
    path "fade_failed_${direction}.tsv", emit: failed
    script:
    def py=(params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    """
    echo '${source_hash}' > code.sha256
    printf '%s\\n' ${fadeJson(expected)} > expected.json
    printf '%s\\n' ${fadeJson(settings)} > settings.json
    ${py} -B ${selectionQuote("${baseDir}/subworkflows/FADE/local/src/fade_merge_tables.py")} \\
        --direction ${direction} --genes ${fadeArgs(genes)} --sites ${fadeArgs(sites)} --audits ${fadeArgs(audits)} \\
        --expected expected.json --settings settings.json
    """
}
