#!/usr/bin/env nextflow
include { FCS_CLOSE } from './fcs_blocks.nf'
include { FCS_CONCAT } from './fcs.nf'

process FCS_RER_OBSERVED {
    label 'process_fcs_batched'
    input:
    tuple val(id), val(size), path(gmts, stageAs: 'gmts/*')
    path stats
    path universe
    output:
    tuple val(id), path('observed.rds'), emit: observed
    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R rer_observed --obs ${stats} --universe ${universe} --gmt-dir gmts --num-g ${params.fcs_min_genes} --max-g ${params.fcs_max_genes ?: 0} --out observed.rds
    """
}
process FCS_RER_NULL_BLOCK {
    label 'process_fcs_batched'
    input:
    tuple val(id), val(size), path(gmts, stageAs: 'gmts/*'), path(observed), path(block)
    path stats
    path universe
    output:
    tuple val(id), path("counts_${id}_${block.baseName}.tsv"), emit: counts
    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R rer_null --matrices ${block} --observed ${observed} --obs ${stats} --universe ${universe} --gmt-dir gmts --num-g ${params.fcs_min_genes} --max-g ${params.fcs_max_genes ?: 0} --out counts_${id}_${block.baseName}.tsv
    """
}
workflow FCS_RER_BLOCKS {
    take:
    stats
    universe
    blocks
    main:
    def sv = stats.first()
    def uv = universe.map { it instanceof List ? it[0] : it }.first()
    def counter = 0
    def batchSize = (params.fcs_batch_size ?: 4) as int
    def batches = Channel.fromPath("${params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"}/*.gmt")
        .collate(batchSize)
        .map { b -> tuple(String.format('fcs_batch_%03d', ++counter), b.size(), b) }
    FCS_RER_OBSERVED(batches, sv, uv)
    def jobs = batches.join(FCS_RER_OBSERVED.out.observed).combine(blocks)
    FCS_RER_NULL_BLOCK(jobs, sv, uv)
    FCS_CLOSE(FCS_RER_OBSERVED.out.observed.join(FCS_RER_NULL_BLOCK.out.counts.groupTuple()))
    FCS_CONCAT(FCS_CLOSE.out.partial.collect())
    emit:
    enrich_file = FCS_CONCAT.out.enrich
}
