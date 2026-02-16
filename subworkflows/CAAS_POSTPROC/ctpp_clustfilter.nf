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
    path "filter_summary.tsv", emit: summary
    path "discarded_summary.tsv", emit: discarded_summary
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    from pathlib import Path
    import re
    
    # Find all filtered files
    filter_files = list(Path('.').glob('*.filtered.*.tsv'))
    
    if not filter_files:
        print("No filtered files found", file=sys.stderr)
        sys.exit(1)
    
    # Check if CAAP mode by examining first file
    first_file = filter_files[0]
    df_sample = pd.read_csv(first_file, sep='\\t', nrows=5)
    has_caap = 'CAAP_Group' in df_sample.columns
    
    print(f"Found {len(filter_files)} filtered files")
    print(f"CAAP mode: {has_caap}")
    
    # Collect per-group counts for each parameter combination
    summary_data = []
    discarded_data = []
    
    for file_path in filter_files:
        # Extract parameters from filename
        fname = file_path.stem
        match = re.search(r'minlen(\\d+)\\.maxcaas(\\d+)', fname)
        if not match:
            continue
        
        minlen = int(match.group(1))
        maxcaas = int(match.group(2))
        param_label = f"{minlen}/{maxcaas}"
        
        # Read file
        df = pd.read_csv(file_path, sep='\\t')
        
        if has_caap and 'CAAP_Group' in df.columns:
            # CAAP mode: count discarded per group
            discarded = df[df['ClusteringFlag'] == 'Discarded']
            
            # Get per-group counts
            group_counts = discarded.groupby('CAAP_Group').size().to_dict()
            
            # Add to summary (one row per parameter with group columns)
            row = {'Parameter': param_label, 'Minlen': minlen, 'Maxcaas': maxcaas}
            row.update(group_counts)
            summary_data.append(row)
            
            # Discarded summary (old format for compatibility)
            groups_str = ','.join(sorted(group_counts.keys()))
            total_count = sum(group_counts.values())
            discarded_data.append({
                'File': fname,
                'DiscardedCount': total_count,
                'Minlen': minlen,
                'Maxcaas': maxcaas,
                'CAAP_Groups': groups_str if groups_str else 'none'
            })
        else:
            # Standard mode: simple count
            count = (df['ClusteringFlag'] == 'Discarded').sum()
            summary_data.append({
                'Parameter': param_label,
                'Minlen': minlen,
                'Maxcaas': maxcaas,
                'DiscardedCount': count
            })
            discarded_data.append({
                'File': fname,
                'DiscardedCount': count,
                'Minlen': minlen,
                'Maxcaas': maxcaas,
                'CAAP_Groups': 'N/A'
            })
    
    # Create DataFrames and save
    if has_caap:
        # For CAAP mode: rows=params, cols=groups
        summary_df = pd.DataFrame(summary_data).fillna(0)
        # Sort by parameters
        summary_df = summary_df.sort_values(['Minlen', 'Maxcaas'])
        # Reorder columns: Parameter, Minlen, Maxcaas, then groups alphabetically
        group_cols = [col for col in summary_df.columns if col not in ['Parameter', 'Minlen', 'Maxcaas']]
        group_cols.sort()
        col_order = ['Parameter', 'Minlen', 'Maxcaas'] + group_cols
        summary_df = summary_df[col_order]
        # Convert group columns to integers
        for col in group_cols:
            summary_df[col] = summary_df[col].astype(int)
    else:
        summary_df = pd.DataFrame(summary_data).sort_values(['Minlen', 'Maxcaas'])
    
    discarded_df = pd.DataFrame(discarded_data).sort_values(['Minlen', 'Maxcaas'])
    
    # Save files
    summary_df.to_csv('filter_summary.tsv', sep='\\t', index=False)
    discarded_df.to_csv('discarded_summary.tsv', sep='\\t', index=False)
    
    print(f"✓ filter_summary.tsv: {len(summary_df)} parameter combinations")
    print(f"✓ discarded_summary.tsv: {len(discarded_df)} entries")
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
