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
    resample_files=\$(find . -type f -name "resample_*.tab" | sort -V)
    
    # Check if we have any files
    if [ -z "\$resample_files" ]; then
        echo "No resample files found" > resample.tab
        exit 0
    fi
    
    # Initialize output with first file (header + data)
    first_file=\$(echo "\$resample_files" | head -n 1)
    cat "\$first_file" > resample.tab
    
    # Append remaining files without their headers
    for file in \$(echo "\$resample_files" | tail -n +2); do
        tail -n +2 "\$file" >> resample.tab
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
    bootstrap_outputs=\$(find . -type f -name "*.output" | sort)
    
    # Check if we have any files
    if [ -z "\$bootstrap_outputs" ]; then
        echo "Gene\tPosition\tCount\tProportion" > bootstrap.tab
        exit 0
    fi
    
    # Initialize output with first file (header + data)
    first_file=\$(echo "\$bootstrap_outputs" | head -n 1)
    cat "\$first_file" > bootstrap.tab
    
    # Append remaining files without their headers
    for file in \$(echo "\$bootstrap_outputs" | tail -n +2); do
        tail -n +2 "\$file" >> bootstrap.tab
    done
    """
}
