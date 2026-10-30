#!/usr/bin/env nextflow
// fcs_blocks.nf — the CAAS FCS on the calibration of SCORING (SCORING_EVIDENCE writes the -log10 p slices).
// Each step is one call of local/src/fcs_blocks.R.
//   FCS_OBSERVED         GMT batch -> observed FCS rows, from the observed -log10 p of all batches
//   FCS_NULL_BLOCK       GMT batch x block x rank -> permulation counts of the block
//   FCS_CLOSE            GMT batch -> counts summed, p.perm, BH and significance (partial table for FCS_CONCAT)
// Batch boundaries of genes come from the universe file in every step: ceil(lines / fcs_block_batch_genes) batches.

nextflow.enable.dsl = 2
include { FCS_CONCAT } from './fcs.nf'

process FCS_OBSERVED {
    tag "${batchID}"
    label 'process_fcs_batched'

    input:
    tuple val(batchID), val(batchSize), path(gmtFiles, stageAs: 'gmts/*')
    path obs_files
    path universe
    path calibration_metadata, stageAs: 'calibration_metadata.json'

    output:
    tuple val(batchID), path("observed.rds"), emit: observed

    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','caas_blocks.R','gene_p.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R observed --obs-dir . --universe ${universe} \\
        --rankings global,top,bottom --gmt-dir gmts --num-g ${params.fcs_min_genes} \\
        --max-g ${params.fcs_max_genes ?: 0} --metadata calibration_metadata.json --out observed.rds
    """
}

process FCS_NULL_BLOCK {
    tag "${batchID} ${block_id} ${rank}"
    label 'process_fcs_batched'

    input:
    tuple val(batchID), val(batchSize), path(gmtFiles, stageAs: 'gmts/*'), path(observed), val(block_id), val(rank), path(gp_files)
    path universe
    val batch_size

    output:
    tuple val(batchID), path("counts_${batchID}_${block_id}_${rank}.tsv"), emit: counts

    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','caas_blocks.R','gene_p.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R null_block --observed ${observed} --slices . \\
        --universe ${universe} --block ${block_id} --batch-size ${batch_size} --rank ${rank} --alt greater \\
        --gmt-dir gmts --num-g ${params.fcs_min_genes} --max-g ${params.fcs_max_genes ?: 0} \\
        --out counts_${batchID}_${block_id}_${rank}.tsv
    """
}

process FCS_CLOSE {
    tag "${batchID}"
    label 'process_fcs_batched'

    input:
    tuple val(batchID), path(observed), path(counts)

    output:
    path "fcs_enrich_partial.tsv", emit: partial

    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','caas_blocks.R','gene_p.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R close --observed ${observed} \\
        --counts \$(ls counts_*.tsv | paste -sd, -) --fdr-thr ${params.fcs_fdr} --pperm-thr ${params.fcs_pperm_thr} \\
        --analytic-levels ${params.fcs_analytic_levels} --out fcs_enrich_partial.tsv
    """
}

// Reuse the single calibration produced by SCORING; no second gene-p computation.
workflow FCS_COMPUTE_CALIBRATED {
    take:
    universe_file
    calibration
    main:
    def ready = calibration.filter { !file("${it}/empty.tsv").exists() }
    // An unavailable calibration (fewer than 2 null cycles, or no genes) has no CAAS FCS. enrich_file then stays empty,
    // so the CAAS FCS report does not run, while AMI and the comparison report still do. The reason is logged here.
    calibration.filter { file("${it}/empty.tsv").exists() }
        .subscribe { d -> log.warn "CAAS FCS skipped: calibration unavailable (fewer than 2 null cycles or no genes), see ${d}/metadata.json" }
    def batchSize = (params.fcs_block_batch_genes ?: 1000) as int
    def gmtDir = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def counter = 0
    def gmtBatchSize = (params.fcs_batch_size ?: 4) as int
    def gmts = Channel.fromPath("${gmtDir}/*.gmt").collate(gmtBatchSize)
        .map { batch -> tuple(String.format('fcs_batch_%03d', ++counter), batch.size(), batch) }
    def obs = ready.map { files("${it}/obs_*.tsv").sort() }
    def gp = ready.map { files("${it}/*_gp_*.rds").sort() }
    def blockIds = ready.flatMap { dir ->
        def meta = new groovy.json.JsonSlurper().parseText(file("${dir}/metadata.json").text)
        meta.blocks instanceof List ? meta.blocks : [meta.blocks]
    }
    def calibration_meta = ready.map { file("${it}/metadata.json") }
    FCS_OBSERVED(gmts, obs, universe_file, calibration_meta)
    def batchObs = gmts.join(FCS_OBSERVED.out.observed.map { id, f -> tuple(id, f) }, by: 0)
    def jobs = batchObs.combine(blockIds).combine(Channel.of('global','top','bottom'))
        .combine(gp.map { [it] }).map { id, size, g, observed, block, rank, fs ->
            tuple(id, size, g, observed, block, rank,
                fs.findAll { it.name.startsWith("${rank}_gp_") && it.name.endsWith("_${block}.rds") })
        }
    FCS_NULL_BLOCK(jobs, universe_file, batchSize)
    FCS_CLOSE(FCS_OBSERVED.out.observed.join(FCS_NULL_BLOCK.out.counts.groupTuple(by:0), by:0)
        .map { id, observed, counts -> tuple(id, observed, counts) })
    // Only a calibration that has the blocks carries metadata.json with the verified gene_stat
    FCS_CONCAT(FCS_CLOSE.out.partial.collect().filter { it }, ready.map { file("${it}/metadata.json") })
    emit:
    enrich_file = FCS_CONCAT.out.enrich
    report_metadata = FCS_CONCAT.out.report_metadata
}
