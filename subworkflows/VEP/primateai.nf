#!/usr/bin/env nextflow

/*
 * PRIMATEAI_MAP
 * ─────────────
 * Maps CAAS variants to PrimateAI-3D pathogenicity scores.
 *
 * Strategy (see map_to_primateai.py for full documentation):
 *   1. Load CAAS file → tag → (caas_pattern, change_side)
 *   2. Load AA2prot output → query → union of changing-side amino acids
 *   3. Load TransVar output → query → (chrom, start, end, transcript, hg38_ref_aa)
 *   4. Build a genomic-position lookup table
 *   5. Stream PrimateAI-3D.hg38.txt.gz and emit rows where:
 *        ref_aa == hg38_ref_aa  AND  alt_aa ∈ changing_aas
 *
 * Output columns:
 *   transvar_query | transvar_transcript | hg38_ref_aa | caas_alt_aas |
 *   [all PrimateAI-3D columns: chr, pos, ref_aa, alt_aa, score_PAI3D,
 *    percentile_PAI3D, refseq, prediction, ...]
 */

process PRIMATEAI_MAP {
    tag "primateai"
    label 'process_long_compute'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/characterization/vep",
               mode: 'copy', overwrite: true,
               pattern: 'primateai_mapped.tsv'

    input:
    path transvar_tsv
    path aa2prot_csv
    path caas_file
    path primateai_db

    output:
    path "primateai_mapped.tsv", emit: primateai_tsv

    stub:
    """
    printf 'transvar_query\ttransvar_transcript\thg38_ref_aa\tcaas_alt_aas\tchr\tpos\tref_aa\talt_aa\tscore_PAI3D\tpercentile_PAI3D\n' > primateai_mapped.tsv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    """
    cp ${local_dir}/map_to_primateai.py .

    if [[ ! -f "${primateai_db}" ]]; then
        echo "WARN Missing PrimateAI database: ${primateai_db}. Skipping PrimateAI mapping." >&2
        touch primateai_mapped.tsv
        exit 0
    fi

    python3 map_to_primateai.py \\
        "${transvar_tsv}" \\
        "${aa2prot_csv}" \\
        "${caas_file}" \\
        "${primateai_db}" \\
        primateai_mapped.tsv
    """
}
