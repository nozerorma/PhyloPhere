#!/usr/bin/env nextflow

/*
#  CT Post-Processing: Cluster filtering (parameter sweep or single filter)
*/

process CAAS_PREPARE_POSTPROC_INPUT {
    tag "prepare_postproc_input"
    publishDir "${params.outdir}/postproc/preprocessed", mode: 'copy', overwrite: true

    input:
    path(discovery_input)

    output:
    path "postproc_discovery_input.tsv", emit: prepared_discovery
    path "removed_patterns_precluster.tsv", emit: removed_patterns

    script:
    def mrca_threshold = params.mrca_posterior_threshold
    """
    python3 - <<'PY'
    import re
    import pandas as pd

    src = "${discovery_input}"
    mrca_threshold = float("${mrca_threshold}")

    def to_bool_series(series: pd.Series) -> pd.Series:
        true_vals = {'true', 't', '1', 'yes', 'y'}
        return (
            series.fillna(False)
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(true_vals)
        )

    # Auto-detect delimiter to support both TSV (legacy discovery) and CSV
    # (ct_disambiguation/caas_convergence_master.csv).
    df = pd.read_csv(src, sep=None, engine='python')

    # Detect disambiguation schema and normalize to postproc schema.
    rename_map = {}
    if 'gene' in df.columns:
        rename_map['gene'] = 'Gene'
    if 'msa_pos' in df.columns:
        rename_map['msa_pos'] = 'Position'
    if 'caap_group' in df.columns:
        rename_map['caap_group'] = 'CAAP_Group'
    if 'pattern_type' in df.columns:
        rename_map['pattern_type'] = 'Pattern'
    if 'pvalue' in df.columns:
        rename_map['pvalue'] = 'Pvalue'

    if rename_map:
        df = df.rename(columns=rename_map)

    if 'Trait' not in df.columns:
        df['Trait'] = 'post_disambiguation'

    if 'CAAP_Group' not in df.columns:
        df['CAAP_Group'] = 'US'

    # Normalize conserved metadata columns for pre-disambiguation pruning.
    if 'is_conserved_meta' not in df.columns and 'IsConserved' in df.columns:
        df['is_conserved_meta'] = df['IsConserved']
    if 'is_conserved_meta' not in df.columns:
        df['is_conserved_meta'] = False
    if 'asr_is_conserved' not in df.columns:
        df['asr_is_conserved'] = False

    df['_is_conserved_meta_bool'] = to_bool_series(df['is_conserved_meta'])
    df['_asr_is_conserved_bool'] = to_bool_series(df['asr_is_conserved'])

    removed_frames = []

    # Keep record of rows removed before clustering for downstream characterization.
    if 'Pattern' in df.columns:
        mask_remove = df['Pattern'].astype(str).str.lower().isin({'ambiguous', 'no_change'})
    else:
        mask_remove = pd.Series([False] * len(df), index=df.index)

    removed_pattern = df.loc[mask_remove].copy()
    if len(removed_pattern) > 0:
        removed_pattern['RemovalReason'] = 'pattern_type_precluster_filter'
        removed_frames.append(removed_pattern)

    cleaned = df.loc[~mask_remove].copy()

    # Remove only conserved-meta rows that failed ASR conserved validation.
    mask_asr_mismatch = cleaned['_is_conserved_meta_bool'] & (~cleaned['_asr_is_conserved_bool'])
    removed_asr = cleaned.loc[mask_asr_mismatch].copy()
    if len(removed_asr) > 0:
        removed_asr['RemovalReason'] = 'conserved_meta_asr_mismatch'
        removed_frames.append(removed_asr)
    cleaned = cleaned.loc[~mask_asr_mismatch].copy()

    # Remove rows with any low-confidence MRCA posterior.
    mrca_cols = [c for c in cleaned.columns if re.fullmatch(r'mrca_\\d+_posterior', str(c))]
    if mrca_cols:
        for col in mrca_cols:
            cleaned[col] = pd.to_numeric(cleaned[col], errors='coerce')

        mask_low_mrca = cleaned[mrca_cols].lt(mrca_threshold).any(axis=1, skipna=True)
        removed_low_mrca = cleaned.loc[mask_low_mrca].copy()
        if len(removed_low_mrca) > 0:
            removed_low_mrca['RemovalReason'] = f'mrca_posterior_below_{mrca_threshold:g}'
            removed_frames.append(removed_low_mrca)
        cleaned = cleaned.loc[~mask_low_mrca].copy()

    # Ensure expected core columns exist for downstream scripts.
    for col in ['Gene', 'Position']:
        if col not in cleaned.columns:
            raise ValueError(f"Required column missing after normalization: {col}")

    cleaned['Position'] = pd.to_numeric(cleaned['Position'], errors='raise').astype(int)

    if removed_frames:
        removed = pd.concat(removed_frames, ignore_index=True)
        if 'Position' in removed.columns:
            removed['Position'] = pd.to_numeric(removed['Position'], errors='coerce')
    else:
        removed = pd.DataFrame(columns=[*df.columns, 'RemovalReason'])

    for helper_col in ['_is_conserved_meta_bool', '_asr_is_conserved_bool']:
        if helper_col in cleaned.columns:
            cleaned = cleaned.drop(columns=[helper_col])
        if helper_col in removed.columns:
            removed = removed.drop(columns=[helper_col])

    cleaned.to_csv('postproc_discovery_input.tsv', sep='\\t', index=False)
    removed.to_csv('removed_patterns_precluster.tsv', sep='\\t', index=False)

    print(f"Input rows: {len(df)}")
    print(f"Removed precluster patterns: {len(removed_pattern)}")
    print(f"Removed conserved-meta ASR mismatch: {len(removed_asr)}")
    if mrca_cols:
        print(f"Removed low MRCA posterior (< {mrca_threshold}): {len(removed_low_mrca)}")
        print(f"MRCA posterior columns used: {', '.join(mrca_cols)}")
    else:
        print("No mrca_*_posterior columns found; MRCA posterior pruning skipped")
    print(f"Rows kept for postproc filtering: {len(cleaned)}")
    PY
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
    python3 ${baseDir}/subworkflows/CT_POSTPROC/local/filter_caas_genes.py \
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
