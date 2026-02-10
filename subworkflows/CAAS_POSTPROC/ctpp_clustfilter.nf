#!/usr/bin/env nextflow

/*
#  CAAS Post-Processing: Cluster filtering (parameter sweep or single filter)
*/

process CAAS_FILTER {
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
    python3 ${baseDir}/subworkflows/CAAS_POSTPROC/local/filter_caas_clusters-param.py \
        -i ${discovery_file} \
        -l ${minlen} \
        -c ${maxcaas} \
        ${params.verbose ? '--verbose' : ''}
    """
}

process CAAS_FILTER_SUMMARY {
    tag "filter_summary"
    publishDir "${params.outdir}/postproc", mode: 'copy', overwrite: true
    
    input:
    path(filter_files)
    
    output:
    path "discarded_summary.tsv", emit: summary
    
    script:
    """
    echo -e "File\tDiscardedCount\tMinlen\tMaxcaas\tCAAP_Groups" > discarded_summary.tsv
    
    for file in *.filtered.*.tsv; do
        # Extract minlen and maxcaas from filename
        # Format: discovery.filtered.minlenX.maxcaasYY.tsv
        base=\$(basename \$file .tsv)
        params=\$(echo \$base | grep -oP 'minlen\\K[0-9]+\\.maxcaas[0-9]+')
        minlen=\$(echo \$params | cut -d'.' -f1)
        maxcaas=\$(echo \$params | cut -d'.' -f2 | sed 's/maxcaas//')
        
        # Check if file has CAAP_Group column
        has_caap_col=\$(head -1 \$file | grep -c "CAAP_Group" || true)
        
        if [ "\$has_caap_col" -gt 0 ]; then
            # CAAP mode: ClusteringFlag is column 4
            # Count discarded entries overall
            count=\$(awk -F'\t' 'NR>1 && \$4=="Discarded" {count++} END {print count+0}' \$file)
            
            # Get unique groups with discarded positions
            groups=\$(awk -F'\t' 'NR>1 && \$4=="Discarded" {print \$3}' \$file | sort -u | tr '\n' ',' | sed 's/,\$//')
            if [ -z "\$groups" ]; then
                groups="none"
            fi
        else
            # CAAS mode: ClusteringFlag is column 3
            count=\$(awk -F'\t' 'NR>1 && \$3=="Discarded" {count++} END {print count+0}' \$file)
            groups="N/A"
        fi
        
        # Append to summary
        echo -e "\${base}\t\${count}\t\${minlen}\t\${maxcaas}\t\${groups}" >> discarded_summary.tsv
    done
    """
}

process CAAS_FILTER_GENES {
    tag "gene_filter:${params.gene_filter_mode}"
    label 'process_low'
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
    
    script:
    def cluster_arg = (params.gene_filter_mode in ['dubious', 'both']) ? "-c ${cluster_file}" : ""
    """
    python3 ${baseDir}/subworkflows/CAAS_POSTPROC/local/filter_caas_genes.py \
        -i ${discovery_file} \
        -l ${gene_ensembl_file} \
        ${cluster_arg} \
        -m ${params.gene_filter_mode} \
        --extreme-percentile ${params.extreme_threshold} \
        --iqr-multiplier ${params.iqr_multiplier} \
        -o filtered_discovery.tsv \
        -s removed_genes_summary.tsv
    """
}

process CAAS_BACKGROUND_CLEANUP {
    tag "bg_cleanup"
    label 'process_low'
    publishDir "${params.outdir}/postproc/cleaned_backgrounds", mode: 'copy', overwrite: true
    
    input:
    path(background_files)
    path(removed_genes_summary)
    
    output:
    path("cleaned_background_*"), emit: cleaned_backgrounds
    
    script:
    """
    python3 ${baseDir}/subworkflows/CAAS_POSTPROC/local/cleanup_background.py \
        -s ${removed_genes_summary} \
        -d . \
        -p "background_genes.*" \
        --extract-trait-from-filename \
        -o .
    """
}
