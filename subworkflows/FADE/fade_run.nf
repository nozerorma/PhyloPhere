#!/usr/bin/env nextflow
// Both execution modes use the same annotation/inference/coverage contract.
include { selectionQuote } from '../SELECTION/selection_contract'
include { fadeJson } from './fade_contract'

process FADE_BATCHED {
    tag "${batchID}"
    label 'process_long_compute'
    cache 'deep'
    publishDir "${params.outdir}/selection/fade/audit", mode:'copy', overwrite:true, pattern:'*.fade.audit.json'

    input:
    tuple val(batchID), val(manifest), path(fastas, stageAs:'fastas/*'), path(trees, stageAs:'trees/*'), path(species, stageAs:'species/*'), val(content_hash)
    val settings
    val source_hash

    output:
    path '*.fade_genes.tsv', emit: genes
    path '*.fade_sites.tsv', emit: sites
    path '*.fade.audit.json', emit: audits
    path '*.hyphy.log', emit: logs, optional:true

    script:
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    def container = (params.use_singularity || params.use_apptainer) ? '--container' : ''
    """
    echo '${source_hash}:${content_hash}' > inputs.sha256
    printf '%s\\n' ${fadeJson(settings)} > settings.json
    cat <<'MANIFEST_EOF' > manifest.tsv
${manifest}MANIFEST_EOF
    ${py} -B ${selectionQuote("${baseDir}/subworkflows/FADE/local/src/run_fade.py")} \\
        --manifest manifest.tsv --settings settings.json --batch-id ${batchID} --workers ${task.cpus} --cpu-per-worker 1 ${container}
    """
}

process FADE_RUN {
    tag "${gene_id}|${direction}"
    label 'process_long_compute'
    cache 'deep'
    publishDir "${params.outdir}/selection/fade/audit", mode:'copy', overwrite:true, pattern:'*.fade.audit.json'

    input:
    tuple val(gene_id), val(direction), val(manifest), path(fastas, stageAs:'fastas/*'), path(trees, stageAs:'trees/*'), path(species, stageAs:'species/*'), val(content_hash)
    val settings
    val source_hash

    output:
    path '*.fade_genes.tsv', emit: genes
    path '*.fade_sites.tsv', emit: sites
    path '*.fade.audit.json', emit: audits
    path '*.hyphy.log', emit: logs, optional:true

    script:
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    def container = (params.use_singularity || params.use_apptainer) ? '--container' : ''
    """
    echo '${source_hash}:${content_hash}' > inputs.sha256
    printf '%s\\n' ${fadeJson(settings)} > settings.json
    cat <<'MANIFEST_EOF' > manifest.tsv
${manifest}MANIFEST_EOF
    ${py} -B ${selectionQuote("${baseDir}/subworkflows/FADE/local/src/run_fade.py")} \\
        --manifest manifest.tsv --settings settings.json --batch-id=${selectionQuote("${gene_id}.${direction}")} --workers 1 --cpu-per-worker ${task.cpus} ${container}
    """
}
