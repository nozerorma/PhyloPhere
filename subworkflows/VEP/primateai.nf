#!/usr/bin/env nextflow

/*
 * PRIMATEAI_MAP
 * ─────────────
 * Maps CAAS variants to PrimateAI-3D pathogenicity scores.
 *
 * Strategy:
 *   1. Load CAAS file → group targets by (Gene, Position)
 *   2. Load MAP files for each Gene → map Position to codon coordinate and infer strand
 *   3. Build a genomic-position lookup table (expand codon to 3 nucleotides)
 *   4. Stream PrimateAI-3D.hg38.txt.gz and emit rows where:
 *        alt_aa ∈ derived_aas  AND  ref_aa matches hg38 reference
 *
 * Output columns:
 *   Gene | Position | hg38_ref_aa | caas_alt_aas | caap_group | scheme_weight |
 *   [all PrimateAI-3D columns: chr, pos, ref_aa, alt_aa, score_PAI3D,
 *    percentile_PAI3D, refseq, prediction, ...]
 */

process PRIMATEAI_MAP {
    tag "primateai"
    label 'process_long_compute'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'primateai_mapped.tsv'

    input:
    path caas_file
    path vep_map_dir
    path primateai_db

    output:
    path "primateai_mapped.tsv", emit: primateai_tsv

    stub:
    """
    printf 'Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaap_group\tscheme_weight\tchr\tpos\tref_aa\talt_aa\tscore_PAI3D\tpercentile_PAI3D\n' > primateai_mapped.tsv
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

    python3 map_to_primateai.py \
        "${caas_file}" \
        "${vep_map_dir}" \
        "${primateai_db}" \
        primateai_mapped.tsv
    """
}
