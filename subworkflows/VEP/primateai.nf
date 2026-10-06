#!/usr/bin/env nextflow
// primateai.nf — Map CAAS positions to PrimateAI-3D pathogenicity scores.
// PhyloPhere | subworkflows/VEP/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  PRIMATEAI_MAP: maps CAAS positions to PrimateAI-3D scores (map_to_primateai.py):
 *    1. group the CAAS rows by (Gene, Position);
 *    2. translate each position to its hg38 codon with the per-gene MAP files
 *       and infer the strand;
 *    3. index the three nucleotides of every codon;
 *    4. stream the PrimateAI-3D table (gzip TSV) and keep the rows whose alt_aa is a
 *       derived residue of the CAAS and (when the ancestral residues are known)
 *       whose ref_aa is an ancestral one.
 *
 *  Consumes:  position_scores.tsv (SCORING), directory of per-gene MAP files,
 *             PrimateAI-3D hg38 table
 *  Produces:  vep/primateai_mapped.tsv, columns Gene, Position, hg38_ref_aa,
 *             caas_alt_aas, caas_change, caap_group, scheme_weight, then every
 *             PrimateAI-3D column (header only when nothing matches; empty when
 *             the database file is missing)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── PrimateAI-3D mapping ─────────────────────────────────────────────────────

process PRIMATEAI_MAP {
    tag "primateai"
    label 'process_long_compute'
    errorStrategy 'ignore'   // a failed annotation never stops the run

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'primateai_mapped.tsv'

    input:
    path position_scores
    path vep_map_dir
    path primateai_db

    output:
    path "primateai_mapped.tsv", emit: primateai_tsv

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    """
    cp ${local_dir}/map_to_primateai.py ${local_dir}/vep_common.py .

    if [[ ! -f "${primateai_db}" ]]; then
        echo "WARN Missing PrimateAI database: ${primateai_db}. Skipping PrimateAI mapping." >&2
        touch primateai_mapped.tsv
        exit 0
    fi

    python3 map_to_primateai.py \
        "${position_scores}" \
        "${vep_map_dir}" \
        "${primateai_db}" \
        primateai_mapped.tsv
    """

    stub:
    """
    printf 'Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaap_group\tscheme_weight\tchr\tpos\tref_aa\talt_aa\tscore_PAI3D\tpercentile_PAI3D\n' > primateai_mapped.tsv
    """
}
