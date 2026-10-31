#!/usr/bin/env nextflow
// Quality masks use every CT candidate, including diagnostic side=none rows.
def ppQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def ppHash(paths) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    paths.each { source ->
        assert source.exists() : "Missing postproc input: ${source}"
        def files = source.isDirectory() ? files("${source}/**/*", type: 'file').sort { it.toString() } : [source]
        files.each { f ->
            md.update((source.isDirectory() ? source.relativize(f).toString() : f.name).getBytes('UTF-8'))
            md.update((byte)0)
            f.toFile().withInputStream { stream ->
                byte[] buffer = new byte[1048576]
                int n
                while ((n = stream.read(buffer)) != -1) md.update(buffer, 0, n)
            }
            md.update((byte)0)
        }
    }
    md.digest().encodeHex().toString()
}

def pp_sources = [
    *files("${baseDir}/subworkflows/CT_POSTPROC/local/**/*.py", type: 'file').findAll { !it.toString().contains('__pycache__') && !it.name.endsWith('.pyc') },
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core/postproc.py"),
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core/scoring_provenance.py"),
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core/columns.py"),
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core/ct_io.py"),
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/utils/io_utils.py"),
    file("${baseDir}/subworkflows/CT/local/modules/design_io.py"),
    file("${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/phylo/species_mapping.py")
]
def pp_source_hash = ppHash(pp_sources.sort { it.toString() })

def clusterParameterGrid(minlens, maxcaases, selected_minlen, selected_maxcaas) {
    (minlens.collectMany { l -> maxcaases.collect { c -> [l,c] } } + [[selected_minlen, selected_maxcaas]]).unique()
}
def clusterFileSuffix(minlen, maxcaas) {
    def token = (new BigDecimal(maxcaas.toString()) * 100).stripTrailingZeros().toPlainString()
    ".filtered.minlen${minlen}.maxcaas${token}.tsv"
}

process CAAS_PREPARE_POSTPROC_INPUT {
    ext source_hash: pp_source_hash
    cache 'deep'
    tag 'prepare_postproc_input'
    publishDir "${params.outdir}/postproc/preprocessed", mode: 'copy', overwrite: true
    input:
    tuple val(observed_hash), path(disambiguation_input, stageAs: 'inputs/master/*'), path(completion, stageAs: 'inputs/completion/*')
    path hyp_pairs, stageAs: 'inputs/design/*'
    path alignment_dir, stageAs: 'inputs/alignment/*'
    val scientific_hash
    path taxid_mapping, stageAs: 'inputs/taxid/*'
    output:
    path 'postproc_disambiguation_input.tsv', emit: prepared_discovery
    path 'removed_patterns_precluster.tsv', emit: removed_patterns
    path 'quality_universe.json', emit: quality_universe
    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def alignment_arg = alignment_dir ? "--alignment ${ppQuote(alignment_dir)} --alignment-format ${ppQuote(params.ali_format ?: 'auto')}" : ''
    def taxid_arg = taxid_mapping ? "--taxid-mapping ${ppQuote(taxid_mapping)}" : ''
    def hyp_arg = hyp_pairs ? "--hyp-pairs ${ppQuote(hyp_pairs)}" : ''
    """
    echo '${task.ext.source_hash}:${observed_hash}:${scientific_hash}' > inputs.sha256
    ${run} ${ppQuote("${baseDir}/subworkflows/CT_POSTPROC/local/src/prepare_postproc_input.py")} \
        --input ${ppQuote(disambiguation_input)} --completion ${ppQuote(completion)} \
        --output postproc_disambiguation_input.tsv --removed-output removed_patterns_precluster.tsv \
        ${alignment_arg} ${hyp_arg} ${taxid_arg}
    """
}

process CT_FILTER {
    ext source_hash: pp_source_hash
    cache 'deep'
    tag "${mode}:${minlen}:${maxcaas}"
    publishDir path: params.caas_postproc_mode == 'exploratory' ? "${params.outdir}/postproc/filter_exploratory/minlen${minlen}_maxcaas${maxcaas}" : "${params.outdir}/postproc/filter_selected", mode: 'copy', overwrite: true
    input:
    tuple val(mode), val(minlen), val(maxcaas), path(discovery_file)
    path map_dir, stageAs: 'inputs/maps/*'
    val map_hash
    output:
    path '*.filtered.*.tsv', emit: filtered_files
    path '*.filtered.*.tsv.scoring.json', emit: filter_contracts
    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def map_arg = map_dir ? "--map-dir ${ppQuote(map_dir)}" : ''
    def count_arg = params.filter_min_positions != null ? "--min-positions ${params.filter_min_positions}" : ''
    """
    echo '${task.ext.source_hash}:${map_hash}' > inputs.sha256
    ${run} ${ppQuote("${baseDir}/subworkflows/CT_POSTPROC/local/src/filter_caas_clusters-param.py")} \
        -i ${ppQuote(discovery_file)} -l ${minlen} -c ${maxcaas} ${count_arg} ${map_arg}
    """
}

process CT_FILTER_SUMMARY {
    ext source_hash: pp_source_hash
    cache 'deep'
    tag 'filter_summary'
    publishDir "${params.outdir}/postproc/summary_statistics", mode: 'copy', overwrite: true, pattern: 'filter_summary.tsv'
    input:
    path filter_files
    output:
    path 'filter_summary.tsv', emit: summary
    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    echo '${task.ext.source_hash}' > source.sha256
    ${run} ${ppQuote("${baseDir}/subworkflows/CT_POSTPROC/local/src/summarize_cluster_filters.py")} --input-dir . --summary-output filter_summary.tsv
    """
}

process CAAS_FILTER_GENES {
    ext source_hash: pp_source_hash
    cache 'deep'
    tag "gene_filter:${params.gene_filter_mode}"
    label 'CT_FILTER'
    publishDir "${params.outdir}/postproc/gene_filtering", mode: 'copy', overwrite: true
    input:
    path discovery_file
    path quality_universe
    path cluster_contract
    path gene_ensembl_file, stageAs: 'inputs/lengths/*'
    path cluster_file
    path map_dir, stageAs: 'inputs/maps/*'
    val map_hash
    output:
    path 'filtered_discovery.tsv', emit: filtered_discovery
    path 'filtered_discovery.tsv.scoring.json', emit: scoring_contract
    path 'removed_genes_summary.tsv', emit: removed_genes
    path 'gene_stats.tsv', emit: gene_stats
    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def map_arg = map_dir ? "--map-dir ${ppQuote(map_dir)}" : ''
    """
    echo '${task.ext.source_hash}:${map_hash}' > inputs.sha256
    ${run} ${ppQuote("${baseDir}/subworkflows/CT_POSTPROC/local/src/filter_caas_genes.py")} \
        -i ${ppQuote(discovery_file)} -l ${ppQuote(gene_ensembl_file)} -c ${ppQuote(cluster_file)} ${map_arg} \
        ${params.remove_caas_clusters ? '--remove-clusters' : ''} -m ${ppQuote(params.gene_filter_mode)} \
        --strata ${ppQuote(params.gene_filter_strata)} --k-far ${ppQuote(params.gene_filter_k_far)} --k-near ${ppQuote(params.gene_filter_k_near)} \
        -o filtered_discovery.tsv -s removed_genes_summary.tsv -g gene_stats.tsv \
        --quality-universe ${ppQuote(quality_universe)} --cluster-contract ${ppQuote(cluster_contract)}
    """
}

process CAAS_BACKGROUND_CLEANUP {
    ext source_hash: pp_source_hash
    cache 'deep'
    tag 'bg_cleanup'
    label 'CT_FILTER'
    publishDir "${params.outdir}/postproc/cleaned_backgrounds", mode: 'copy', overwrite: true
    input:
    path global_background_file, stageAs: 'inputs/background/*'
    path removed_genes_summary
    output:
    path 'cleaned_background_main.txt', emit: cleaned_background_main
    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    echo '${task.ext.source_hash}' > source.sha256
    ${run} ${ppQuote("${baseDir}/subworkflows/CT_POSTPROC/local/src/cleanup_background.py")} \
        -s ${ppQuote(removed_genes_summary)} -g ${ppQuote(global_background_file)} -o .
    """
}
