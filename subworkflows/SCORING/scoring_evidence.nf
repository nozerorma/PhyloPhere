// Shared calibration, independent of GMT/FCS/AMI switches.
process SCORING_EVIDENCE {
    label 'process_fcs_batched'
    publishDir path: "${params.outdir}/scoring", mode: 'copy', overwrite: true, pattern: '*.tsv'
    publishDir path: "${params.outdir}/scoring", mode: 'copy', overwrite: true, pattern: '*_lists'
    publishDir path: "${params.outdir}/scoring", mode: 'copy', overwrite: true, pattern: 'gene_evidence_metadata.json'
    input:
    path genes, stageAs: 'base/gene_scores.tsv'
    path positions, stageAs: 'base/position_scores.tsv'
    path stats, stageAs: 'base/fcs_stats.tsv'
    path universe, stageAs: 'universe.txt'
    path live_blocks, stageAs: 'live/*'
    path perms, stageAs: 'null_input'
    output:
    path 'gene_scores.tsv', emit: gene_scores
    path 'position_scores.tsv', emit: position_scores
    path 'fcs_stats.tsv', emit: fcs_stats
    path 'gene_evidence.tsv', emit: gene_evidence
    path 'gene_evidence_metadata.json', emit: metadata
    path 'gene_lists', emit: gene_lists
    path 'position_lists', emit: position_lists
    path 'calibration', emit: calibration
    script:
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def source_files = ['SCORING/local/src/scoring_evidence.R','SCORING/local/src/caas_block_contract.R','ENRICHMENT/local/src/fcs_prepare_blocks.R','ENRICHMENT/local/src/gene_p.R','ENRICHMENT/local/src/caas_blocks.R']
    def hash = java.security.MessageDigest.getInstance('SHA-256').digest(source_files.collect { file("${baseDir}/subworkflows/${it}").text }.join('\n').getBytes('UTF-8')).encodeHex().toString()
    """
    # evidence sources: ${hash}
    ${rs} ${baseDir}/subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R caas ${params.gene_stat} ${params.caas_cycle_block_size ?: 500} ${params.caas_permulation_enrichment == false ? 'disabled' : 'enabled'}
    ${rs} ${baseDir}/subworkflows/SCORING/local/src/scoring_evidence.R --genes ${genes} --positions ${positions} --stats ${stats} --universe ${universe} --blocks blocks --batch ${params.fcs_block_batch_genes ?: 1000} --gene-stat ${params.gene_stat} --nominal ${params.evidence_p_threshold ?: 0.05} --fdr ${params.evidence_fdr_threshold ?: 0.05}
    """
}
