#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_MATRIX {
    tag "RER_MATRIX|${params.traitname}"
    label 'process_medium'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_objects", mode: 'copy'

    input:
    path trait_file
    path gene_trees_file

    val source_hash
    output:
    path "${params.traitname}.RERmatrix.output", emit: matrix

    script:
    def input_hash = rerHash([trait_file, gene_trees_file])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_matrix.R")} ${rerQuote(trait_file)} ${rerQuote(gene_trees_file)} ${rerQuote("${params.traitname}.RERmatrix.output")} ${rerQuote(params.rer_minsp)}
    """
}
