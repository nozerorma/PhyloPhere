#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_TABLES {
    tag "RER_TABLES|${params.traitname}"
    label 'process_low'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy', pattern: '{rerconverge_summary_*.tsv,rer_run_metadata.json,rer_inputs}'

    input:
    path result_rds
    path null_input

    val source_hash
    output:
    path 'validated_rer.output', emit: validated
    path "rerconverge_summary_${params.traitname}.tsv", emit: summary_tsv
    path 'rer_run_metadata.json', emit: metadata
    path 'rer_inputs', emit: inputs

    script:
    def input_hash = rerHash([result_rds, null_input])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_tables.R")} ${rerQuote(result_rds)} ${rerQuote(null_input)} ${rerQuote(params.traitname)} ${rerQuote(params.evidence_p_threshold == null ? 0.05 : params.evidence_p_threshold)} ${rerQuote(params.evidence_fdr_threshold == null ? 0.05 : params.evidence_fdr_threshold)}
    """
}
