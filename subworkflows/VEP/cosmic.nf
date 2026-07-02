#!/usr/bin/env nextflow

/*
 * COSMIC_MAP
 * ──────────
 * Maps CAAS variants to COSMIC Mutant Census somatic mutations.
 */

process COSMIC_MAP {
    tag "cosmic"
    label 'process_long_compute'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'cosmic_scores.tsv'

    input:
    path caas_file
    path vep_map_dir
    path cosmic_db

    output:
    path "cosmic_scores.tsv", emit: cosmic_tsv

    stub:
    """
    printf 'Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaas_change\tcaap_group\tscheme_weight\tCHROMOSOME\tGENOME_START\tGENOMIC_WT_ALLELE\tGENOMIC_MUT_ALLELE\tMUTATION_AA\tMUTATION_DESCRIPTION\tMUTATION_SOMATIC_STATUS\n' > cosmic_scores.tsv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    """
    cp ${local_dir}/map_to_cosmic.py .

    if [[ ! -f "${cosmic_db}" ]]; then
        echo "WARN Missing COSMIC database: ${cosmic_db}. Skipping COSMIC mapping." >&2
        touch cosmic_scores.tsv
        exit 0
    fi

    python3 map_to_cosmic.py \
        "${caas_file}" \
        "${vep_map_dir}" \
        "${cosmic_db}" \
        cosmic_scores.tsv
    """
}
