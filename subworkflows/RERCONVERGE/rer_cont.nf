#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_CONT {
    tag "RER_CONT|${params.traitname}"
    label 'process_medium'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_results", mode: 'copy'

    input:
    path trait_file
    path rer_master_tree
    path rer_matrix

    val source_hash
    output:
    path "${params.traitname}.char2path.output", emit: char2path
    path "${params.traitname}.continuous.output", emit: continuous_output
    path "${params.traitname}.continuous.perms.rds", emit: perms_output, optional: true

    script:
    def input_hash = rerHash([trait_file, rer_master_tree, rer_matrix])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def batches = params.rer_perm_batches == null ? 0 : params.rer_perm_batches
    def per_batch = params.rer_perms_per_batch == null ? 100 : params.rer_perms_per_batch

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    rm -rf -- ${rerQuote("${params.traitname}.continuous.perms.rds")} ${rerQuote("${params.traitname}.continuous.perms.rds.gene_slices")}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/continuous_rer.R")} ${rerQuote(trait_file)} ${rerQuote(rer_master_tree)} ${rerQuote("${params.traitname}.char2path.output")} ${rerQuote(rer_matrix)} ${rerQuote("${params.traitname}.continuous.output")} ${rerQuote(params.rer_minsp)} ${rerQuote(params.winsorizeRER)} ${rerQuote(params.winsorizeTrait)} ${rerQuote(batches)} ${rerQuote(per_batch)} cc ${rerQuote(params.rer_transform ?: 'auto')} ${rerQuote(params.rer_batch_seed)}
    """
}
