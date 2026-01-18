/*
 * Concatenation processes for CAAS tools outputs
 * Combines multiple output files, keeping only the first header
 */

process CONCAT_DISCOVERY {
    tag "Concatenating discovery outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(discovery_files)

    output:
    path("discovery.tab"), emit: discovery_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Find all .output files and sort them
    mapfile -t discovery_files < <(find . -type f -name "*.output" | sort)
    
    # Check if we have any files
    if [ \${#discovery_files[@]} -eq 0 ]; then
        echo "Gene\tPosition\tAncestral\tDerived\tType\tForeground\tBackground" > discovery.tab
        exit 0
    fi
    
    # Copy first file completely (header + data)
    cat "\${discovery_files[0]}" > discovery.tab
    
    # Append remaining files without their headers
    for ((i=1; i<\${#discovery_files[@]}; i++)); do
        tail -n +2 "\${discovery_files[\$i]}" >> discovery.tab
    done
    """
}

process CONCAT_RESAMPLE {
    tag "Concatenating resample outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(resample_dir)

    output:
    path("resample.tab"), emit: resample_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Find all resample_*.tab files and sort them numerically
    mapfile -t resample_files < <(find . -type f -name "resample_*.tab" | sort -V)
    
    # Check if we have any files
    if [ \${#resample_files[@]} -eq 0 ]; then
        echo "No resample files found" > resample.tab
        exit 0
    fi
    
    # Copy first file completely (header + data)
    cat "\${resample_files[0]}" > resample.tab
    
    # Append remaining files without their headers
    for ((i=1; i<\${#resample_files[@]}; i++)); do
        tail -n +2 "\${resample_files[\$i]}" >> resample.tab
    done
    """
}

process CONCAT_BOOTSTRAP {
    tag "Concatenating bootstrap outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(bootstrap_files)

    output:
    path("bootstrap.tab"), emit: bootstrap_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    # Find all .output files and sort them
    mapfile -t bootstrap_files < <(find . -type f -name "*.output" | sort)
    
    # Check if we have any files
    if [ \${#bootstrap_files[@]} -eq 0 ]; then
        echo "Gene\tPosition\tCount\tProportion" > bootstrap.tab
        exit 0
    fi
    
    # Copy first file completely (header + data)
    cat "\${bootstrap_files[0]}" > bootstrap.tab
    
    # Append remaining files without their headers
    for ((i=1; i<\${#bootstrap_files[@]}; i++)); do
        tail -n +2 "\${bootstrap_files[\$i]}" >> bootstrap.tab
    done
    """
}
