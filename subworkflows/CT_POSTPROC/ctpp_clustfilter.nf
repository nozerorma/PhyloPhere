#!/usr/bin/env nextflow

/*
#  CT Post-Processing: Cluster filtering (parameter sweep or single filter)
*/

process CAAS_PREPARE_POSTPROC_INPUT {
    tag "prepare_postproc_input"
    publishDir "${params.outdir}/postproc/preprocessed", mode: 'copy', overwrite: true

    input:
    path(disambiguation_input)

    output:
    path "postproc_disambiguation_input.tsv", emit: prepared_discovery
    path "removed_patterns_precluster.tsv", emit: removed_patterns

    script:
    def mrca_threshold = params.ct_disambig_posterior_threshold
<<<<<<< HEAD
<<<<<<< HEAD
=======
    def use_all_mrca_filter = params.use_all_mrca_filter ? "true" : "false"
>>>>>>> asr
=======
>>>>>>> convergence
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/prepare_postproc_input.py \
        --input ${disambiguation_input} \
        --mrca-threshold ${mrca_threshold} \
        --output postproc_disambiguation_input.tsv \
        --removed-output removed_patterns_precluster.tsv
    """
}

process CT_FILTER {
    tag "${mode}:${minlen}x${maxcaas_int}"
    publishDir(
        path: params.caas_postproc_mode == 'exploratory' ? 
            "${params.outdir}/postproc/filter_${mode}/minlen${minlen}_maxcaas${maxcaas_int}" :
            "${params.outdir}/postproc/filter_selected",
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(mode), val(minlen), val(maxcaas), path(discovery_file)
    
    output:
    path "*.filtered.*.tsv", emit: filtered_files
    
    script:
    maxcaas_int = (maxcaas * 100).toInteger()
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/filter_caas_clusters-param.py \
        -i ${discovery_file} \
        -l ${minlen} \
        -c ${maxcaas} \
        ${params.verbose ? '--verbose' : ''}
    """
}

process CT_FILTER_SUMMARY {
    tag "filter_summary"
    publishDir "${params.outdir}/postproc", mode: 'copy', overwrite: true
    
    input:
    path(filter_files)
    
    output:
    path "filter_summary.tsv", emit: summary
    path "discarded_summary.tsv", emit: discarded_summary

    script:
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/summarize_cluster_filters.py \
        --input-dir . \
        --summary-output filter_summary.tsv \
        --discarded-output discarded_summary.tsv
    """
}

process CAAS_FILTER_GENES {
    tag "gene_filter:${params.gene_filter_mode}"
    label 'CT_FILTER'
    publishDir "${params.outdir}/postproc/gene_filtering", mode: 'copy', overwrite: true
    
    when:
    params.gene_filter_mode != 'none'
    
    input:
    path(discovery_file)
    path(gene_ensembl_file)
    path(cluster_file)
    
    output:
    path "filtered_discovery.tsv", emit: filtered_discovery
    path "removed_genes_summary.tsv", emit: removed_genes
    path "gene_stats.tsv", emit: gene_stats, optional: true
    
    script:
    def cluster_arg = (params.gene_filter_mode in ['dubious', 'both']) ? "-c ${cluster_file}" : ""
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/filter_caas_genes.py \
        -i ${discovery_file} \
        -l ${gene_ensembl_file} \
        ${cluster_arg} \
        -m ${params.gene_filter_mode} \
        --extreme-percentile ${params.extreme_threshold} \
        --iqr-multiplier ${params.iqr_multiplier} \
        -o filtered_discovery.tsv \
        -s removed_genes_summary.tsv \
        -g gene_stats.tsv
    """
}

process CAAS_BACKGROUND_CLEANUP {
    tag "bg_cleanup"
    label 'CT_FILTER'
    publishDir "${params.outdir}/postproc/cleaned_backgrounds", mode: 'copy', overwrite: true
    
    input:
    path(global_background_file)
    path(removed_genes_summary)
    
    output:
    path("cleaned_background_*"), emit: cleaned_backgrounds
    path("cleaned_background_main.txt"), emit: cleaned_background_main
    
    script:
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/cleanup_background.py \
        -s ${removed_genes_summary} \
        -g ${global_background_file} \
        -o .
    """
}
