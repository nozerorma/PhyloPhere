/*
 * Concatenation processes for CAAS tools outputs
 * Combines multiple output files, keeping only the first header
 */

process CONCAT_DISCOVERY {
    tag "Concatenating discovery outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    val(discovery_dir)

    output:
    path("discovery.tab"), emit: discovery_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_DISCOVERY DEBUG ==="
    echo "Working directory: \$(pwd)"
    echo "Discovery directory: ${discovery_dir}"
    echo ""
    
    # Find all .output files in the published directory and sort them
    mapfile -t discovery_files < <(find "${discovery_dir}" -type f -name "*.output" | sort)
    
    echo "Found \${#discovery_files[@]} discovery files:"
    printf '%s\n' "\${discovery_files[@]}"
    echo ""
    
    # Check if we have any files
    if [ \${#discovery_files[@]} -eq 0 ]; then
        echo "WARNING: No discovery files found - creating header-only file"
        echo "Gene\tPosition\tAncestral\tDerived\tType\tForeground\tBackground" > discovery.tab
        exit 0
    fi
    
    echo "First file preview (\${discovery_files[0]}):"
    head -5 "\${discovery_files[0]}" || echo "ERROR: Cannot read first file"
    echo "Line count: \$(wc -l < "\${discovery_files[0]}")"
    echo ""
    
    # Copy first file completely (header + data)
    cat "\${discovery_files[0]}" > discovery.tab
    
    # Append remaining files without their headers
    for ((i=1; i<\${#discovery_files[@]}; i++)); do
        echo "Appending file \$((i+1))/\${#discovery_files[@]}: \${discovery_files[\$i]} (\$(wc -l < "\${discovery_files[\$i]}") lines)"
        tail -n +2 "\${discovery_files[\$i]}" >> discovery.tab
    done
    
    echo ""
    echo "Final concatenated file line count: \$(wc -l < discovery.tab)"
    echo "Final file preview:"
    head -10 discovery.tab
    """
}

process CONCAT_RESAMPLE {
    tag "Concatenating resample outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    val(resample_dir)

    output:
    path("resample.tab"), emit: resample_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_RESAMPLE DEBUG ==="
    echo "Working directory: \$(pwd)"
    echo "Resample directory: ${resample_dir}"
    echo ""
    
    # Find all resample_*.tab files in the published directory and sort them numerically
    mapfile -t resample_files < <(find "${resample_dir}" -type f -name "resample_*.tab" | sort -V)
    
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

process CONCAT_BOOTSTRAP {
    tag "Concatenating bootstrap outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    val(bootstrap_dir)

    output:
    path("bootstrap.tab"), emit: bootstrap_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_BOOTSTRAP DEBUG ==="
    echo "Working directory: \$(pwd)"
    echo "Bootstrap directory: ${bootstrap_dir}"
    echo ""
    
    # Find all .output files in the published directory and sort them
    mapfile -t bootstrap_files < <(find "${bootstrap_dir}" -type f -name "*.bootstraped.output" | sort)
    
    echo "Found \${#bootstrap_files[@]} bootstrap files:"
    printf '%s\n' "\${bootstrap_files[@]}"
    echo ""
    
    # Check if we have any files
    if [ \${#bootstrap_files[@]} -eq 0 ]; then
        echo "WARNING: No bootstrap files found - creating header-only file"
        echo "Gene\tPosition\tCount\tProportion" > bootstrap.tab
        exit 0
    fi
    
    echo "First file preview (\${bootstrap_files[0]}):"
    head -5 "\${bootstrap_files[0]}" || echo "ERROR: Cannot read first file"
    echo "Line count: \$(wc -l < "\${bootstrap_files[0]}")"
    echo ""
    
    # Copy first file completely
    cat "\${bootstrap_files[0]}" > bootstrap.tab
    
    # Append all remaining files (bootstrap files don't have headers, all lines are data)
    for ((i=1; i<\${#bootstrap_files[@]}; i++)); do
        echo "Appending file \$((i+1))/\${#bootstrap_files[@]}: \${bootstrap_files[\$i]} (\$(wc -l < "\${bootstrap_files[\$i]}") lines)"
        cat "\${bootstrap_files[\$i]}" >> bootstrap.tab
    done
    
    echo ""
    echo "Final concatenated file line count: \$(wc -l < bootstrap.tab)"
    echo "Final file preview:"
    head -10 bootstrap.tab
    """
}
