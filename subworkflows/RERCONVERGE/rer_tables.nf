#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_TABLES {
    tag "RER_TABLES|${params.traitname}"
    label 'process_low'
    errorStrategy 'terminate'
    cache true
    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy', pattern: '{rerconverge_summary_*.tsv,rer_run_metadata.json,rer_gene_identifiers.tsv,rer_inputs,normalized_null}'

    input:
    path result_rds
    path null_input
    path gene_reference, stageAs: 'gene_reference/*'

    val source_hash
    output:
    path 'validated_rer.output', emit: validated
    path "rerconverge_summary_${params.traitname}.tsv", emit: summary_tsv
    path 'rer_run_metadata.json', emit: metadata
    path 'rer_inputs', emit: inputs
    path 'normalized_null', emit: normalized_null, optional: true

    script:
    def input_hash = rerHash([result_rds, null_input, gene_reference])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_tables.R")} ${rerQuote(result_rds)} ${rerQuote(null_input)} ${rerQuote(params.traitname)} ${rerQuote(params.evidence_p_threshold == null ? 0.05 : params.evidence_p_threshold)} ${rerQuote(params.evidence_fdr_threshold == null ? 0.05 : params.evidence_fdr_threshold)} ${rerQuote(gene_reference)} ${rerQuote("${baseDir}/subworkflows/SCORING/local/src/scoring_io.R")}
    """
}
