#!/usr/bin/env nextflow
// fcs_blocks.nf — the CAAS FCS on blocks of the null (WP6). Each step is one call of local/src/fcs_blocks.R.
//   FCS_BLOCK_TRANSPOSE  block -> one slice per batch of genes and rank (global, top, bottom)
//   FCS_GENE_P           rank x batch of genes -> -log10 p slices of every block, observed -log10 p of the batch
//   FCS_OBSERVED         GMT batch -> observed FCS rows, from the observed -log10 p of all batches
//   FCS_NULL_BLOCK       GMT batch x block x rank -> permulation counts of the block
//   FCS_CLOSE            GMT batch -> counts summed, p.perm, BH and significance (partial table for FCS_CONCAT)
// Batch boundaries of genes come from the universe file in every step: ceil(lines / fcs_block_batch_genes) batches.

nextflow.enable.dsl = 2
include { FCS_CONCAT } from './fcs.nf'

process FCS_BLOCK_TRANSPOSE {
    tag "${block_id}"
    label 'process_fcs_batched'

    input:
    tuple val(block_id), path(matrices)
    path universe
    val batch_size

    output:
    path "{global,top,bottom}_b*.rds", emit: slices   // the slices only, not the staged matrices.rds

    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','caas_blocks.R','gene_p.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R transpose --matrices ${matrices} \\
        --universe ${universe} --batch ${batch_size} --rank all --block ${block_id} --gene-stat ${params.gene_stat} --out .
    """
}

process FCS_GENE_P {
    tag "${rank} batch ${k}"
    label 'process_fcs_batched'

    input:
    tuple val(rank), val(k), val(obs_col), path(slices)
    path stats_file

    output:
    path "*_gp_*.rds", emit: gp
    path "obs_*.tsv",  emit: obs

    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_hash = java.security.MessageDigest.getInstance('SHA-256').digest(
        ['fcs_blocks.R','caas_blocks.R','gene_p.R','fcs_enrich.R'].collect {
            file("${baseDir}/subworkflows/ENRICHMENT/local/src/${it}").text
        }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # FCS sources: ${source_hash}
    blocks=\$(ls ${rank}_b\$(printf '%04d' ${k})_*.rds | sed -E "s/^${rank}_b[0-9]{4}_//; s/[.]rds\$//" | paste -sd, -)
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R gene_p --batch ${k} --rank ${rank} \\
        --blocks "\$blocks" --slices . --obs ${stats_file} --obs-col ${obs_col} --out .
    """
}

process FCS_OBSERVED {
    tag "${batchID}"
    label 'process_fcs_batched'

    input:
    tuple val(batchID), val(batchSize), path(gmtFiles, stageAs: 'gmts/*')
    path obs_files
    path universe

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
        --max-g ${params.fcs_max_genes ?: 0} --out observed.rds
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
    nb=\$(( (\$(grep -c . ${universe}) + ${batch_size} - 1) / ${batch_size} ))
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_blocks.R null_block --observed ${observed} --slices . \\
        --universe ${universe} --block ${block_id} --batches "\$(seq -s, 1 \$nb)" --rank ${rank} --alt greater \\
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

workflow FCS_COMPUTE_BLOCKS {
    take:
    stats_file        // path: gene + score_global, score_top, score_bottom (selected gene_stat)
    universe_file     // path: gene universe
    blocks            // channel of tuple(block_id, matrices.rds)

    main:
    def batchSize = (params.fcs_block_batch_genes ?: 1000) as int
    def gmt_dir_resolved = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def gmtBatchSize = (params.fcs_batch_size ?: 4) as int
    def counter = 0
    def gmtBatches = Channel.fromPath("${gmt_dir_resolved}/*.gmt")
        .collate(gmtBatchSize)
        .map { batch -> tuple(String.format('fcs_batch_%03d', ++counter), batch.size(), batch) }

    def universe_v = universe_file.first()
    def stats_v    = stats_file.first()
    def nb_ch = universe_v.flatMap { f -> (1..((f.readLines().findAll { it.trim() }.size() + batchSize - 1).intdiv(batchSize))) }

    FCS_BLOCK_TRANSPOSE(blocks, universe_v, batchSize)
    def slices_v = FCS_BLOCK_TRANSPOSE.out.slices.collect()

    def rank_cols = Channel.of(['global', 'score_global'], ['top', 'score_top'], ['bottom', 'score_bottom'])
    def gp_jobs = rank_cols.combine(nb_ch).map { rank, col, k -> tuple(rank, k, col) }
        .combine(slices_v.map { files -> [files] }).map { rank, k, col, sl -> tuple(rank, k, col, sl.findAll { it.name.startsWith(String.format('%s_b%04d_', rank, k)) }) }   // the list stays one item
    FCS_GENE_P(gp_jobs, stats_v)
    def gp_v  = FCS_GENE_P.out.gp.collect()
    def obs_v = FCS_GENE_P.out.obs.collect()

    FCS_OBSERVED(gmtBatches, obs_v.filter { it }.first(), universe_v)

    // every GMT batch with its observed rows, crossed with each block and rank
    def batch_obs = gmtBatches.join(FCS_OBSERVED.out.observed.map { id, f -> tuple(id, f) }, by: 0)
        .map { id, size, gmts, obs -> tuple(id, size, gmts, obs) }
    def nul_jobs = batch_obs.combine(blocks.map { it[0] }).combine(Channel.of('global', 'top', 'bottom'))
        .combine(gp_v.map { files -> [files] })
        .map { id, size, gmts, obs, block, rank, files ->
            tuple(id, size, gmts, obs, block, rank, files.findAll { it.name.startsWith("${rank}_gp_") && it.name.endsWith("_${block}.rds") }) }
    FCS_NULL_BLOCK(nul_jobs, universe_v, batchSize)

    def counts_by_batch = FCS_NULL_BLOCK.out.counts.groupTuple(by: 0)
    FCS_CLOSE(FCS_OBSERVED.out.observed.join(counts_by_batch, by: 0).map { id, obs, cts -> tuple(id, obs, cts) })

    FCS_CONCAT(FCS_CLOSE.out.partial.collect().filter { it })

    emit:
    enrich_file = FCS_CONCAT.out.enrich
    partials = FCS_CLOSE.out.partial.collect()   // concatenate with FCS_CONCAT (fcs.nf)
}

// Reuse the single calibration produced by SCORING; no second gene-p computation.
workflow FCS_COMPUTE_CALIBRATED {
    take:
    universe_file
    calibration
    main:
    def ready = calibration.filter { !file("${it}/empty.tsv").exists() }
    def empty = calibration.filter { file("${it}/empty.tsv").exists() }.map { file("${it}/empty.tsv") }
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
    FCS_OBSERVED(gmts, obs, universe_file)
    def batchObs = gmts.join(FCS_OBSERVED.out.observed.map { id, f -> tuple(id, f) }, by: 0)
    def jobs = batchObs.combine(blockIds).combine(Channel.of('global','top','bottom'))
        .combine(gp.map { [it] }).map { id, size, g, observed, block, rank, fs ->
            tuple(id, size, g, observed, block, rank,
                fs.findAll { it.name.startsWith("${rank}_gp_") && it.name.endsWith("_${block}.rds") })
        }
    FCS_NULL_BLOCK(jobs, universe_file, batchSize)
    FCS_CLOSE(FCS_OBSERVED.out.observed.join(FCS_NULL_BLOCK.out.counts.groupTuple(by:0), by:0)
        .map { id, observed, counts -> tuple(id, observed, counts) })
    FCS_CONCAT(FCS_CLOSE.out.partial.collect().filter { it })
    emit:
    enrich_file = FCS_CONCAT.out.enrich.mix(empty)
}
