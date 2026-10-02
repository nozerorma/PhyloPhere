/*
 * Concatenation of the resample output of the CT tools
 */

process CONCAT_RESAMPLE {
    tag "Concatenating resample outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(resample_dir)

    output:
    path("resample.tab"), emit: resample_concat
    // The per-cycle audit trail written by permulations.R (tier, pair count,
    // overall Dunn, and the permulated trait values the FG/BG were selected on).
    // Published alongside resample.tab so the pool can be checked from the
    // results directory rather than only from the Nextflow work dir. Optional:
    // resample_from / precomputed runs stage a directory that has no manifest.
    path("permulation_manifest.tsv"), emit: resample_manifest, optional: true

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Carry the manifest through when the staged directory has one.
    if [ -f "${resample_dir}/permulation_manifest.tsv" ]; then
        cp "${resample_dir}/permulation_manifest.tsv" permulation_manifest.tsv
        echo "Manifest carried through: \$(wc -l < permulation_manifest.tsv) lines"
    else
        echo "No permulation_manifest.tsv in staged dir (precomputed or legacy resample)"
    fi

    echo "=== CONCAT_RESAMPLE ==="
    echo "Working directory: \$(pwd)"
    echo "Staged directory: ${resample_dir}"
    echo "Contents of staged directory:"
    ls -la ${resample_dir}/
    echo ""
    
    # Find all resample_*.tab files in the staged directory and sort them numerically
    mapfile -t resample_files < <(find ${resample_dir}/ -type f -name "resample_*.tab" | sort -V)
    
    echo "Found \${#resample_files[@]} resample files:"
    printf '%s\n' "\${resample_files[@]}"
    echo ""
    
    # Check if we have any files
    if [ \${#resample_files[@]} -eq 0 ]; then
        echo "WARNING: No resample files found - creating placeholder file"
        echo "No resample files found" > resample.tab
        exit 0
    fi
    
    echo "First file preview (\${resample_files[0]}):"
    head -5 "\${resample_files[0]}" || echo "ERROR: Cannot read first file"
    echo "Line count: \$(wc -l < "\${resample_files[0]}")"
    echo ""
    
    # Copy first file completely
    cat "\${resample_files[0]}" > resample.tab
    
    # Append all remaining files (resample files don't have headers, all lines are data)
    for ((i=1; i<\${#resample_files[@]}; i++)); do
        echo "Appending file \$((i+1))/\${#resample_files[@]}: \${resample_files[\$i]} (\$(wc -l < "\${resample_files[\$i]}") lines)"
        cat "\${resample_files[\$i]}" >> resample.tab
    done
    
    echo ""
    echo "Final concatenated file line count: \$(wc -l < resample.tab)"
    echo "Final file preview:"
    head -10 resample.tab
    """
}
