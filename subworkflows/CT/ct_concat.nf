// ct_concat.nf — Concatenate the partitioned resample output into a single resample.tab.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CONCAT_RESAMPLE: joins the resample_*.tab files that RESAMPLE writes (version-sorted
 *  by name, no header line in them) into resample.tab for reporting, and carries
 *  permulation_manifest.tsv and permulation_harvest.tsv through when the staged directory has them.
 *
 *  Consumes:  resample directory (RESAMPLE output, or a directory given by resample_from)
 *  Produces:  resample.tab, permulation_manifest.tsv and permulation_harvest.tsv (optional), published to
 *             caastools/
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Resample concatenation ───────────────────────────────────────────────────

process CONCAT_RESAMPLE {
    tag "Concatenating resample outputs"
    publishDir "${params.outdir}/caastools", mode: 'copy'

    input:
    path(resample_dir)

    output:
    path("resample.tab"), emit: resample_concat
    // Per-cycle audit trail written by permulations.R (tier, pair count, overall Dunn, and the
    // permulated trait values the FG/BG were selected on). It is published next to resample.tab so
    // the pool can be checked from the results directory. Optional: a directory given by
    // resample_from may have no manifest.
    path("permulation_manifest.tsv"), emit: resample_manifest, optional: true
    // Audit of the harvest (draws, rejections by reason, acceptance by PSS tolerance, pool and capacity of a
    // sample of the draws), written by permulations.R; the "Null harvest" tab of the scoring report reads it.
    path("permulation_harvest.tsv"), emit: resample_harvest, optional: true

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

    if [ -f "${resample_dir}/permulation_harvest.tsv" ]; then
        cp "${resample_dir}/permulation_harvest.tsv" permulation_harvest.tsv
    fi

    echo "=== CONCAT_RESAMPLE ==="
    echo "Working directory: \$(pwd)"
    echo "Staged directory: ${resample_dir}"
    echo "Contents of staged directory:"
    ls -la ${resample_dir}/
    echo ""
    
    # Version-sorted, so that resample_010.tab follows resample_009.tab
    mapfile -t resample_files < <(find ${resample_dir}/ -type f -name "resample_*.tab" | sort -V)
    
    echo "Found \${#resample_files[@]} resample files:"
    printf '%s\n' "\${resample_files[@]}"
    echo ""
    
    # Without resample files a placeholder is written and the process ends successfully
    if [ \${#resample_files[@]} -eq 0 ]; then
        echo "WARNING: No resample files found - creating placeholder file"
        echo "No resample files found" > resample.tab
        exit 0
    fi
    
    echo "First file preview (\${resample_files[0]}):"
    head -5 "\${resample_files[0]}" || echo "ERROR: Cannot read first file"
    echo "Line count: \$(wc -l < "\${resample_files[0]}")"
    echo ""
    
    # The first file starts resample.tab
    cat "\${resample_files[0]}" > resample.tab
    
    # The resample files have no header line, so the others are appended whole
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
