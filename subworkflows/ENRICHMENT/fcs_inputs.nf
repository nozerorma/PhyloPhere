#!/usr/bin/env nextflow
include { FCS_COMPUTE_BLOCKS } from './fcs_blocks.nf'
include { FCS_RER_BLOCKS } from './fcs_rer_blocks.nf'

// Resolve live block directories, a saved block directory, or the compatibility RDS.
process FCS_PREPARE_BLOCKS {
    label 'process_fcs_batched'
    input:
    path live_blocks, stageAs: 'live/*'
    path perms, stageAs: 'null_input'
    val source
    output:
    path 'blocks/*.rds', emit: blocks, optional: true
    path 'empty.tsv', emit: empty, optional: true
    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def code = file("${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R").text
    def hash = java.security.MessageDigest.getInstance('SHA-256').digest(code.getBytes('UTF-8')).encodeHex().toString()
    """
    # resolver: ${hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R ${source} ${params.gene_stat ?: 'size_adj_max'} ${source == 'caas' ? (params.caas_cycle_block_size ?: 500) : (params.rer_perms_per_batch ?: 100)} ${source == 'caas' && params.caas_permulation_enrichment == false ? 'disabled' : 'enabled'}
    """
}

workflow FCS_CAAS {
    take:
    stats
    universe
    perms
    live_blocks
    main:
    FCS_PREPARE_BLOCKS(live_blocks.ifEmpty { [] }, perms, 'caas')
    def blocks = FCS_PREPARE_BLOCKS.out.blocks.flatten().map { f -> tuple(f.baseName, f) }
    FCS_COMPUTE_BLOCKS(stats, universe, blocks)
    emit:
    enrich_file = FCS_COMPUTE_BLOCKS.out.enrich_file.mix(FCS_PREPARE_BLOCKS.out.empty)
}

workflow FCS_RER {
    take:
    stats
    universe
    perms
    main:
    FCS_PREPARE_BLOCKS(Channel.value([]), perms, 'rer')
    FCS_RER_BLOCKS(stats, universe, FCS_PREPARE_BLOCKS.out.blocks.flatten())
    emit:
    enrich_file = FCS_RER_BLOCKS.out.enrich_file.mix(FCS_PREPARE_BLOCKS.out.empty)
}
