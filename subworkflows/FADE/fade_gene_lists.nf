#!/usr/bin/env nextflow
include { selectionQuote; selectionHash } from '../SELECTION/selection_contract'
include { fadeSources } from './fade_contract'
process FADE_GENE_LISTS_TASK {
    tag "${direction}"
    label 'process_low'
    cache 'deep'
    publishDir "${params.outdir}/selection/fade/${direction}/gene_lists", mode:'copy', overwrite:true
    input:
    val direction
    tuple path(summary_tsv), val(input_hash)
    val source_hash
    output:
    val direction, emit: direction
    path '*.txt', emit: gene_lists
    script:
    def py=(params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    """
    echo '${source_hash}:${input_hash}' > inputs.sha256
    ${py} -B ${selectionQuote("${baseDir}/subworkflows/FADE/local/src/fade_gene_lists.py")} \\
        --direction ${direction} --summary ${selectionQuote(summary_tsv)}
    """
}

workflow FADE_GENE_LISTS {
    take:
    direction
    summary_ch
    main:
    def inputs=summary_ch.map { f -> tuple(f,selectionHash([f])) }
    def result=FADE_GENE_LISTS_TASK(direction,inputs,fadeSources('lists'))
    emit:
    direction=result.direction
    gene_lists=result.gene_lists
}
