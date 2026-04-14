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
    def transvarCfg = params.vep_transvar_cfg ?: ''
    def transvarReference = params.vep_transvar_reference ?: ''
    """
    # Extract unique GENE:p.N query strings from column 4 of the AA2prot table
    cut -f4 "${aa2prot_csv}" | sort -u > positions_input.txt

    staged_db_dir="\$(cd "${transvar_db}" && pwd)"

    # Resolve TransVar config with explicit precedence:
    #  1) user-provided --vep_transvar_cfg (if it exists)
    #  2) staged DB transvar.config
    #  3) staged DB transvar.cfg
    # If found, stage as a task-local temporary config and export TRANSVAR_CFG.
    cfg_src=""
    if [[ -n "${transvarCfg}" ]]; then
        if [[ -f "${transvarCfg}" ]]; then
            cfg_src="${transvarCfg}"
        else
            echo "WARN Requested TransVar config not found: ${transvarCfg}. Trying staged defaults." >&2
        fi
    fi

    if [[ -z "\${cfg_src}" && -f "\${staged_db_dir}/transvar.config" ]]; then
        cfg_src="\${staged_db_dir}/transvar.config"
    fi
    if [[ -z "\${cfg_src}" && -f "\${staged_db_dir}/transvar.cfg" ]]; then
        cfg_src="\${staged_db_dir}/transvar.cfg"
    fi

    if [[ -n "\${cfg_src}" ]]; then
        cp "\${cfg_src}" transvar.cfg
        export TRANSVAR_CFG="\$(pwd)/transvar.cfg"
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

    # TransVar does not accept --dbdir on the panno subcommand. Prefer an
    # explicit reference FASTA when provided, then attempt to resolve it from
    # 'transvar config --refversion'. If unresolved, still attempt panno and
    # let TransVar use TRANSVAR_CFG/default behavior.
    ref_path="${transvarReference}"
    if [[ -n "\${ref_path}" && ! -f "\${ref_path}" ]]; then
        echo "WARN vep_transvar_reference points to a missing file: \${ref_path}. Ignoring it." >&2
        ref_path=""
    fi

    if [[ -z "\${ref_path}" ]]; then
        ref_path=\$("\${tv_bin}" config --refversion "${refver}" 2>/dev/null | awk -F': ' '
            tolower(\$1) ~ /reference/ { print \$2; exit }
        ')
    fi

    ref_args=()
    if [[ -n "\${ref_path}" ]]; then
        ref_args=(--reference "\${ref_path}")
    else
        echo "WARN TransVar could not resolve a reference FASTA for refversion ${refver}. Continuing without --reference and relying on TRANSVAR_CFG/default config." >&2
    fi

    ccds_db="\${staged_db_dir}/${refver}.ccds.txt.transvardb"
    if [[ ! -f "\${ccds_db}" ]]; then
        echo "WARN TransVar CCDS database not found: \${ccds_db}. Skipping TransVar annotation." >&2
        touch transvar.tsv
        exit 0
    fi

    export TRANSVAR_DOWNLOAD_DIR="\${staged_db_dir}"

    if ! "\${tv_bin}" panno \\
        -l positions_input.txt \\
        --ccds "\${ccds_db}" \\
        --refversion "${refver}" \\
        "\${ref_args[@]}" \\
        --noheader \\
        --longest \\
        > transvar.tsv; then
        echo "WARN TransVar panno failed for refversion ${refver}. Skipping TransVar annotation." >&2
        touch transvar.tsv
        exit 0
    fi

    # Guarantee the output file exists even if TransVar produced nothing
    [[ -s transvar.tsv ]] || touch transvar.tsv
    """
}
