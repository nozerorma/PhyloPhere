#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_GENE_LISTS {
    tag "RER_GENE_LISTS|${params.traitname}"
    label 'process_low'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/gene_lists", mode: 'copy'

    input:
    path summary_tsv
    path gene_scores

    val source_hash
    output:
    path '*.txt', emit: gene_lists
    path 'fcs_stats.tsv', emit: fcs_stats

    script:
    def input_hash = rerHash([summary_tsv, gene_scores])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_gene_lists.R")} ${rerQuote(summary_tsv)} ${rerQuote(gene_scores)} ${rerQuote(params.evidence_p_threshold == null ? 0.05 : params.evidence_p_threshold)}
    """
}
