#!/usr/bin/env nextflow
include { selectionQuote; selectionHash } from '../SELECTION/selection_contract'
include { fadeSources } from './fade_contract'
process FADE_RENDER {
    tag "${direction}"
    label 'process_reporting'
    cache 'deep'
    publishDir "${params.outdir}/selection/fade/${direction}", mode:'copy', overwrite:true
    publishDir "${params.outdir}/html_reports", mode:'copy', overwrite:true, pattern:'*.html'
    input:
    val direction
    tuple path(gene_counts_file), path(sites_file), path(fg_list_file), path(metadata_file), val(input_hash)
    val source_hash
    output:
    path "6.FADE_report_${direction}.html", emit: report
    path 'fade_inputs', emit: saved_inputs
    script:
    def r=(params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    echo '${source_hash}:${input_hash}' > inputs.sha256
    cp ${selectionQuote("${baseDir}/subworkflows/FADE/local/6.FADE_report.Rmd")} .
    cp ${selectionQuote("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")} .
    ${r} ${selectionQuote("${baseDir}/subworkflows/FADE/local/src/render_fade.R")} \\
        ${direction} ${selectionQuote(gene_counts_file)} ${selectionQuote(sites_file)} ${selectionQuote(fg_list_file)} ${selectionQuote(metadata_file)}
    """
}

workflow FADE_REPORT {
    take:
    direction
    counts
    sites
    foreground
    metadata
    main:
    def inputs=counts.combine(sites).combine(foreground).combine(metadata).map { c,s,f,m ->
        tuple(c,s,f,m,selectionHash([c,s,f,m].findAll { !it.name.startsWith('NO_') }))
    }
    def result=FADE_RENDER(direction,inputs,fadeSources('report'))
    emit:
    report=result.report
    saved_inputs=result.saved_inputs
}
