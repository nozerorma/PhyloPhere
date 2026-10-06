#!/usr/bin/env nextflow
// cosmic.nf — Map CAAS positions to COSMIC Mutant Census somatic mutations.
// PhyloPhere | subworkflows/VEP/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  COSMIC_MAP: converts each CAAS position to its hg38 codon through the per-gene
 *  MAP files and keeps the COSMIC missense mutations of that codon whose
 *  ancestral→derived amino-acid change matches the CAAS (map_to_cosmic.py).
 *
 *  Consumes:  position_scores.tsv (SCORING), directory of per-gene MAP files,
 *             COSMIC Mutant Census GRCh38 table (gzip TSV)
 *  Produces:  vep/cosmic_scores.tsv (header only when nothing matches; empty when
 *             the database file is missing)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── COSMIC mapping ───────────────────────────────────────────────────────────

process COSMIC_MAP {
    tag "cosmic"
    label 'process_long_compute'
    errorStrategy 'ignore'   // a failed annotation never stops the run

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'cosmic_scores.tsv'

    input:
    path position_scores
    path vep_map_dir
    path cosmic_db

    output:
    path "cosmic_scores.tsv", emit: cosmic_tsv

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    """
    cp ${local_dir}/map_to_cosmic.py ${local_dir}/vep_common.py .

    if [[ ! -f "${cosmic_db}" ]]; then
        echo "WARN Missing COSMIC database: ${cosmic_db}. Skipping COSMIC mapping." >&2
        touch cosmic_scores.tsv
        exit 0
    fi

    python3 map_to_cosmic.py \
        "${position_scores}" \
        "${vep_map_dir}" \
        "${cosmic_db}" \
        cosmic_scores.tsv
    """

    stub:
    """
    printf 'Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaas_change\tcaap_group\tscheme_weight\tCHROMOSOME\tGENOME_START\tGENOMIC_WT_ALLELE\tGENOMIC_MUT_ALLELE\tMUTATION_AA\tMUTATION_DESCRIPTION\tMUTATION_SOMATIC_STATUS\n' > cosmic_scores.tsv
    """
}
