/*
 * Concatenation processes for CT tools outputs
 * Combines multiple output files, keeping only the first header
 */

process CONCAT_DISCOVERY {
    tag "Concatenating discovery outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(discovery_files, stageAs: "discovery_*")

    output:
    path("discovery.tab"), emit: discovery_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_DISCOVERY ==="
    echo "Working directory: \$(pwd)"
    echo "Staged files:"
    ls -la
    echo ""
    
    # Find all discovery_* files in the staged directory (files or symlinks)
    mapfile -t discovery_files < <(find . -maxdepth 1 -name "discovery_*" ! -name ".*" | sort)
    
    echo "Found \${#discovery_files[@]} discovery files:"
    printf '%s\n' "\${discovery_files[@]}"
    echo ""

    # Check if we have any files
    if [ \${#discovery_files[@]} -eq 0 ]; then
        echo "WARNING: No discovery files found - creating empty discovery table with expected schema"
        # Keep schema aligned with CT_SIGNIFICATION expectations (Pattern + Pvalue required)
        echo "Gene\tMode\tCAAP_Group\tTrait\tPosition\tSubstitution\tEncoded\tPvalue\tPattern\tFFGN\tFBGN\tGFG\tGBG\tMFG\tMBG\tFFG\tFBG\tMS\tConservedPair\tConservedPairs" > discovery.tab
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

process CONCAT_BACKGROUND {
    tag "Concatenating background outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(background_files, stageAs: "background_*")

    output:
    path("background.output"), emit: background_concat
    path("background_genes.output"), emit: background_genes

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_BACKGROUND ==="
    echo "Working directory: \$(pwd)"
    echo "Staged files:"
    ls -la
    echo ""
    
    # Find all background_* files in the staged directory (files or symlinks)
    mapfile -t background_files < <(find . -maxdepth 1 -name "background_*" ! -name ".*" | sort)
    
    echo "Found \${#background_files[@]} background files:"
    printf '%s\n' "\${background_files[@]}"
    echo ""

    # Check if we have any files
    if [ \${#background_files[@]} -eq 0 ]; then
        echo "WARNING: No background files found - creating header-only file"
        echo "Gene\tPosition" > background.output
        touch background_genes.output
        exit 0
    fi

    echo "First file preview (\${background_files[0]}):"
    head -5 "\${background_files[0]}" || echo "ERROR: Cannot read first file"
    echo "Line count: \$(wc -l < "\${background_files[0]}")"
    echo ""

    # Concatenate all background files (no headers to strip)
    cat "\${background_files[@]}" > background.output

    echo "Final background file line count: \$(wc -l < background.output)"
    echo "Final file preview:"
    head -10 background.output
    
    # Generate background_genes.output: unique gene list where Position is not empty
    echo ""
    echo "Generating background_genes.output..."
    awk -F'\t' 'NF >= 2 && \$2 != "" && \$2 != "Position" {print \$1}' background.output | sort -u > background_genes.output
    echo "Unique genes with positions: \$(wc -l < background_genes.output)"
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

process CONCAT_BOOTSTRAP {
    tag "Concatenating bootstrap outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(bootstrap_files, stageAs: "bootstrap_*")

    output:
    path("bootstrap.tab"), emit: bootstrap_concat

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "=== CONCAT_BOOTSTRAP ==="
    echo "Working directory: \$(pwd)"
    echo "Staged files:"
    ls -la
    echo ""
    
    # Find all bootstrap_* files in the staged directory (files or symlinks)
    mapfile -t bootstrap_files < <(find . -maxdepth 1 -name "bootstrap_*" ! -name ".*" | sort)
    
    echo "Found \${#bootstrap_files[@]} bootstrap files:"
    printf '%s\n' "\${bootstrap_files[@]}"
    echo ""
    
    # Check if we have any files
    if [ \${#bootstrap_files[@]} -eq 0 ]; then
        echo "WARNING: No bootstrap files found - creating header-only file"
        echo "Gene@Position\tCAAP_Group\tCount\tTotal\tProportion" > bootstrap.tab
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
