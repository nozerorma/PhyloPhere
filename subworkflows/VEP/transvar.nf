#!/usr/bin/env nextflow

/*
 * TRANSVAR_ANNO
 * ─────────────
 * Runs TransVar panno in batch on the GENE:p.N protein-position query strings
 * produced by PROT2AA.
 *
 * Flags used:
 *   --ccds       use CCDS canonical transcripts only
 *   --longest    return one annotation per query (canonical/longest CDS)
 *   --noheader   skip the header line in the output
 *   --dbdir      path to the TransVar index directory (dat/transvar/)
 *
 * Output columns (tab-separated, no header):
 *   input | transcript | gene | strand | coordinates(gDNA/cDNA/protein) | region | info
 *
 * transvar must be available on PATH.  When running inside the phylophere
 * container or conda env it is pre-installed.
 */

process TRANSVAR_ANNO {
    tag "transvar"
    label 'process_vep'

    publishDir path: "${params.outdir}/characterization/vep",
               mode: 'copy', overwrite: true,
               pattern: 'transvar.tsv'

    input:
    path aa2prot_csv
    path transvar_db

    output:
    path "transvar.tsv", emit: transvar_tsv

    stub:
    """
    printf 'input\ttranscript\tgene\tstrand\tcoordinates(gDNA/cDNA/protein)\tregion\tinfo\n' > transvar.tsv
    """

    script:
    def refver = params.vep_refversion ?: 'hg38'
    """
    # Extract unique GENE:p.N query strings from column 4 of the AA2prot table
    cut -f4 "${aa2prot_csv}" | sort -u > positions_input.txt

    transvar panno \\
        -l positions_input.txt \\
        --ccds \\
        --refversion "${refver}" \\
        --noheader \\
        --longest \\
        --dbdir "${transvar_db}" \\
        > transvar.tsv 2>/dev/null || true

    # Guarantee the output file exists even if TransVar produced nothing
    [[ -s transvar.tsv ]] || touch transvar.tsv
    """
}
