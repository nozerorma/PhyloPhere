#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_REPORT {
    tag "rer_report|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy', pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', pattern: '*.html'
    input:
    path continuous_output
    val source_hash
    output:
    path '5.RERconverge_report.html', emit: report
    script:
    def local_dir = "${baseDir}/subworkflows/RERCONVERGE/local"
    def style = file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")
    def input_hash = rerHash([continuous_output])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def top = params.rer_top_n_labels == null ? 15 : params.rer_top_n_labels
    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    cp ${rerQuote(local_dir+'/5.RERconverge_report.Rmd')} .
    cp ${rerQuote(local_dir+'/rer_contract.R')} .
    cp ${rerQuote(style)} report_style.R
    ${rs} ${rerQuote(local_dir+'/render_rer.R')} ${rerQuote(continuous_output)} ${rerQuote(params.traitname)} ${rerQuote(top)}
    """
}
