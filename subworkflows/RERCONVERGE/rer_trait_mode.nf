#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_TRAIT_MODE {
    tag "RER_TRAIT_MODE|${params.traitname}"
    label 'process_low'
    errorStrategy 'terminate'
    cache 'deep'

    input:
    path trait_file

    val source_hash
    output:
    path 'rer_mode.txt', emit: mode

    script:
    def input_hash = rerHash([trait_file])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_trait_mode.R")} ${rerQuote(trait_file)} ${rerQuote(params.rer_trait_mode ?: 'auto')} rer_mode.txt
    """
}
