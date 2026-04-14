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
    errorStrategy 'ignore'

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
    def transvarBin = params.vep_transvar_bin ?: ''
    def transvarCfg = params.vep_transvar_cfg ?: '${projectDir}/subworkflows/VEP/dat/transvar/transvar.cfg'
    """
    # Extract unique GENE:p.N query strings from column 4 of the AA2prot table
    cut -f4 "${aa2prot_csv}" | sort -u > positions_input.txt

    if [[ -n "${transvarCfg}" ]]; then
        export TRANSVAR_CFG="${transvarCfg}"
    fi

    if [[ -x "${transvarBin}" ]]; then
        tv_bin="${transvarBin}"
    elif command -v transvar >/dev/null 2>&1; then
        tv_bin="\$(command -v transvar)"
    else
        echo "WARN TransVar executable not found (checked: ${transvarBin} and PATH). Skipping TransVar annotation." >&2
        touch transvar.tsv
        exit 0
    fi

    # TransVar does not accept --dbdir on the panno subcommand. Resolve the
    # active reference FASTA, then pass the staged CCDS .transvardb path
    # explicitly on the command line so we do not inherit stale paths from any
    # existing user/site TransVar configuration.
    ref_path=\$("\${tv_bin}" config --refversion "${refver}" 2>/dev/null | awk -F': ' '/^Reference:/ {print \$2; exit}')
    if [[ -z "\${ref_path}" ]]; then
        echo "WARN TransVar could not resolve a reference FASTA for refversion ${refver}. Skipping TransVar annotation." >&2
        touch transvar.tsv
        exit 0
    fi

    staged_db_dir="\$(cd "${transvar_db}" && pwd)"
    ccds_db="\${staged_db_dir}/${refver}.ccds.txt.transvardb"
    if [[ ! -f "\${ccds_db}" ]]; then
        echo "WARN TransVar CCDS database not found: \${ccds_db}. Skipping TransVar annotation." >&2
        touch transvar.tsv
        exit 0
    fi

    export TRANSVAR_DOWNLOAD_DIR="\${staged_db_dir}"

    "\${tv_bin}" panno \\
        -l positions_input.txt \\
        --ccds "\${ccds_db}" \\
        --refversion "${refver}" \\
        --reference "\${ref_path}" \\
        --noheader \\
        --longest \\
        > transvar.tsv

    # Guarantee the output file exists even if TransVar produced nothing
    [[ -s transvar.tsv ]] || touch transvar.tsv
    """
}
