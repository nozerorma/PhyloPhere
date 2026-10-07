#!/usr/bin/env nextflow
// ctpp_clustfilter.nf — Input preparation, cluster filter and gene filter of CT post-processing.
// PhyloPhere | subworkflows/CT_POSTPROC/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_POSTPROC filters: prepare the disambiguation table, flag clustered positions for
 *  one (minlen, maxcaas) pair or for a sweep of pairs, summarize the sweep, remove
 *  extreme and dubious genes, and clean the background gene list. Called from
 *  workflows/ct_postproc.nf.
 *
 *  Consumes:  disambiguation master table, optional alignments and contrast design,
 *             gene annotation (lengths), global background gene list
 *  Produces:  postproc_disambiguation_input.tsv, per-pair cluster files
 *             (*.filtered.minlen*.maxcaas*.tsv), filter_summary.tsv, discarded_summary.tsv,
 *             filtered_discovery.tsv, removed_genes_summary.tsv, gene_stats.tsv and
 *             cleaned_background_*.txt, under postproc/
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Helper functions ─────────────────────────────────────────────────────────

// (minlen, maxcaas) pairs of the exploratory sweep: minlen_values x maxcaas_values, plus the selected pair
// (filter_minlen, filter_maxcaas) when the grid does not contain it, so the gene filter always has its cluster file.
def clusterParameterGrid(minlens, maxcaases, selected_minlen, selected_maxcaas) {
    def combos = []
    minlens.each { l -> maxcaases.each { c -> combos << [l, c] } }
    if (!combos.any { it[0] == selected_minlen && it[1] == selected_maxcaas }) {
        combos << [selected_minlen, selected_maxcaas]
    }
    return combos
}

// Suffix of the cluster file that filter_caas_clusters-param.py writes for a (minlen, maxcaas) pair.
def clusterFileSuffix(minlen, maxcaas) {
    return ".filtered.minlen${minlen}.maxcaas${(maxcaas * 100).toInteger()}.tsv"
}


// ── Processes ────────────────────────────────────────────────────────────────

// Normalizes the disambiguation master table (prepare_postproc_input.py).
process CAAS_PREPARE_POSTPROC_INPUT {
    tag "prepare_postproc_input"
    publishDir "${params.outdir}/postproc/preprocessed", mode: 'copy', overwrite: true

    input:
    path(disambiguation_input)
    path(hyp_pairs)   // contrast_hypotheses_pairs.tsv, or the NO_HYP_PAIRS sentinel

    output:
    path "postproc_disambiguation_input.tsv", emit: prepared_discovery
    path "removed_patterns_precluster.tsv", emit: removed_patterns

    script:
    // Optional extant-species residue tally: it needs the alignment directory (params.alignment)
    // and the contrast design (contrast_hypotheses_pairs.tsv: the top / bottom species of each
    // hypothesis), so it does not depend on FADE. Without them prepare_postproc_input.py leaves
    // the tally columns empty.
    def ali_dir   = params.alignment ?: ''
    def ali_fmt   = params.ali_format ?: 'fasta'
    def ali_flag  = ali_dir ? "--alignment '${ali_dir}' --alignment-format '${ali_fmt}'" : ''
    def hyp_flag  = (hyp_pairs.name =~ /^NO_/) ? '' : "--hyp-pairs '${hyp_pairs}'"
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/src/prepare_postproc_input.py \
        --input ${disambiguation_input} \
        --output postproc_disambiguation_input.tsv \
        --removed-output removed_patterns_precluster.tsv \
        ${ali_flag} ${hyp_flag}
    """
}

// Flags the clustered positions of one (minlen, maxcaas) pair (filter_caas_clusters-param.py).
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
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/src/filter_caas_clusters-param.py \
        -i ${discovery_file} \
        -l ${minlen} \
        -c ${maxcaas} \
        ${params.caas_map_dir ? "--map-dir '${params.caas_map_dir}'" : ''} \
        --verbose
    """
}

// Counts the discarded positions of every cluster file (summarize_cluster_filters.py).
process CT_FILTER_SUMMARY {
    tag "filter_summary"
    publishDir "${params.outdir}/postproc/summary_statistics", mode: 'copy', overwrite: true, pattern: "filter_summary.tsv"
    publishDir "${params.outdir}/postproc/gene_filtering", mode: 'copy', overwrite: true, pattern: "discarded_summary.tsv"

    input:
    path(filter_files)

    output:
    path "filter_summary.tsv", emit: summary
    path "discarded_summary.tsv", emit: discarded_summary

    script:
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/src/summarize_cluster_filters.py \
        --input-dir . \
        --summary-output filter_summary.tsv \
        --discarded-output discarded_summary.tsv
    """
}

// Removes extreme and dubious genes, and the flagged positions when params.remove_caas_clusters is set (filter_caas_genes.py).
process CAAS_FILTER_GENES {
    tag "gene_filter:${params.gene_filter_mode}"
    label 'CT_FILTER'
    publishDir "${params.outdir}/postproc/gene_filtering", mode: 'copy', overwrite: true

    input:
    path(discovery_file)
    path(gene_ensembl_file)
    path(cluster_file)

    output:
    path "filtered_discovery.tsv", emit: filtered_discovery
    path "removed_genes_summary.tsv", emit: removed_genes
    path "gene_stats.tsv", emit: gene_stats, optional: true

    script:
    def cluster_arg = cluster_file ? "-c ${cluster_file}" : ""
    def remove_clusters_flag = params.remove_caas_clusters ? "--remove-clusters" : ""
    """
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/src/filter_caas_genes.py \
        -i ${discovery_file} \
        -l ${gene_ensembl_file} \
        ${cluster_arg} \
        ${remove_clusters_flag} \
        -m ${params.gene_filter_mode} \
        --extreme-percentile ${params.extreme_threshold} \
        --iqr-multiplier ${params.iqr_multiplier} \
        -o filtered_discovery.tsv \
        -s removed_genes_summary.tsv \
        -g gene_stats.tsv
    """
}

// Drops the removed genes from the background gene list (cleanup_background.py).
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
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/src/cleanup_background.py \
        -s ${removed_genes_summary} \
        -g ${global_background_file} \
        -o .
    """
}
