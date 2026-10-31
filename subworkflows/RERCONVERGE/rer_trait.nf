#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_TRAIT {
    tag "RER_TRAIT|${params.traitname}"
    label 'process_low'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_traits", mode: 'copy'

    input:
    path my_traitfile

    val source_hash
    output:
    path "${params.traitname}.polished.output", emit: polished
    path "${params.traitname}.trait_type.output", emit: trait_type

    script:
    def input_hash = rerHash([my_traitfile])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/build_rer_trait.R")} ${rerQuote(my_traitfile)} ${rerQuote(params.sp_colname)} ${rerQuote(params.traitname)} ${rerQuote("${params.traitname}.polished.output")} ${rerQuote("${params.traitname}.trait_type.output")} ${rerQuote(params.n_trait ?: '')} ${rerQuote(params.c_trait ?: '')} ${rerQuote(params.rer_transform ?: 'auto')} ${rerQuote(params.rer_trait_mode ?: 'auto')}
    """
}
