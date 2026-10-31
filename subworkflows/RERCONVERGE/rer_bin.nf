#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_BIN {
    tag "RER_BIN|${params.traitname}"
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
    path "${params.traitname}.fg_paths.output", emit: fg_paths
    path "${params.traitname}.binary.output", emit: binary_output
    path "${params.traitname}.binary.perms.rds", emit: perms_output, optional: true

    script:
    def input_hash = rerHash([trait_file, rer_master_tree, rer_matrix])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def batches = params.rer_perm_batches == null ? 0 : params.rer_perm_batches
    def per_batch = params.rer_perms_per_batch == null ? 100 : params.rer_perms_per_batch

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    rm -rf -- ${rerQuote("${params.traitname}.binary.perms.rds")} ${rerQuote("${params.traitname}.binary.perms.rds.gene_slices")}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/binary_rer.R")} ${rerQuote(trait_file)} ${rerQuote(rer_master_tree)} ${rerQuote("${params.traitname}.fg_paths.output")} ${rerQuote(rer_matrix)} ${rerQuote("${params.traitname}.binary.output")} ${rerQuote(params.rer_minsp)} ${rerQuote(params.rer_min_pos == null ? 2 : params.rer_min_pos)} ${rerQuote(params.winsorizeRER)} ${rerQuote(params.rer_binary_clade ?: 'all')} ${rerQuote(batches)} ${rerQuote(per_batch)} ${rerQuote(params.rer_batch_seed)}
    """
}
