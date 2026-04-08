#!/usr/bin/env nextflow

/*
 * PROT2COORD
 * ──────────
 * Translates CAAS protein positions to genomic coordinates using the human
 * CDS alignment, GFF annotation, and per-gene codon-filtering tracking files.
 *
 * For each gene in the CAAS file the script:
 *   1. Reads the selected/removed codon indices from the TRACK HTML file.
 *   2. Maps the human protein position to nucleotide position via the CDS alignment.
 *   3. Uses the GFF to resolve the genomic coordinate (strand-aware).
 *
 * Output columns: Gene | Position | tag | chr | current_coord | real_coord
 */

process PROT2COORD {
    tag "prot2coord"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/characterization/vep",
               mode: 'copy', overwrite: true,
               pattern: 'aa2nuc_global.csv'

    input:
    path caas_file
    path cds_dir
    path track_dir
    path hs_cds_gz
    path gff_file
    path gene_equiv

    output:
    path "aa2nuc_global.csv", emit: aa2nuc_csv

    stub:
    """
    printf 'Gene\tPosition\ttag\tchr\tcurrent_coord\treal_coord\n' > aa2nuc_global.csv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    """
    cp ${local_dir}/Extract_genomic_coordinates.py .

    mkdir -p out
    failed=0

    genes=\$(awk -F'\\t' '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (tolower(\$i) == "gene") { col=i; break }
            }
            next
        }
        col && \$col != "" { print \$col }
    ' "${caas_file}" | sort -u)

    while IFS= read -r gene; do
        [[ -z "\$gene" ]] && continue
        track=\$(find -L "${track_dir}" \\( -name "\${gene}.html" -o -name "\${gene}.*.html" \\) -print -quit 2>/dev/null || true)
        cds=\$(find -L "${cds_dir}" \\( \
            -name "\${gene}.Homo*.fasta" -o -name "\${gene}.*.fasta" -o \
            -name "\${gene}.Homo*.fa"    -o -name "\${gene}.*.fa"    -o \
            -name "\${gene}.Homo*.fna"   -o -name "\${gene}.*.fna" \
        \\) -print -quit 2>/dev/null || true)
        if [[ -n "\$track" && -n "\$cds" ]]; then
            python3 Extract_genomic_coordinates.py \\
                "\$cds" "${hs_cds_gz}" "\$track" "${gff_file}" \\
                "${caas_file}" "out/\${gene}_out.csv" "${gene_equiv}" "\$gene" \\
                || { echo "ERROR prot2coord: script failed for \${gene}" >&2; failed=1; }
        else
            [[ -z "\$track" ]] && echo "SKIP prot2coord \${gene}: tracking file not found" >&2
            [[ -z "\$cds" ]]   && echo "SKIP prot2coord \${gene}: CDS fasta not found" >&2
        fi
    done <<< "\$genes"

    (( failed == 0 )) || exit 1

    shopt -s nullglob
    out_files=(out/*.csv)
    if (( \${#out_files[@]} )); then
        cat "\${out_files[@]}" > aa2nuc_global.csv
    else
        touch aa2nuc_global.csv
    fi
    shopt -u nullglob
    """
}
