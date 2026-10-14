#!/usr/bin/env nextflow
// ensembl_vep.nf — Annotate CAAS amino-acid changes with Ensembl VEP consequences.
// PhyloPhere | subworkflows/VEP/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ENSEMBL_VEP_ANNOTATE: annotates the ancestral→derived amino-acid change of every
 *  CAAS position with the consequence prediction of the Ensembl variant_effect_predictor
 *  (the `vep` command line tool, not this pipeline's --vep toggle). It does not
 *  depend on PrimateAI-3D or COSMIC.
 *
 *  Steps: build_vep_hgvs.py writes one protein-level HGVS identifier per change of every
 *  (gene, position, side) (the ancestral and derived residues are the ASR ones of
 *  position_scores.tsv and differ between the two sides, so no reference proteome is needed),
 *  `vep` annotates them, and join_vep_output.py attaches the result to Gene, Position and
 *  side. The output has one row per transcript; the reports reduce it to one
 *  consequence per (Gene, Position, side). The command requests --canonical and, for vep 105
 *  or later, --mane, so that the CANONICAL, MANE and MANE_SELECT columns flag the canonical
 *  and MANE Select transcripts that reduction prefers. The input is protein-level, so the
 *  consequences are those of a coding substitution (missense, start lost, splice region,
 *  stop retained and the like); UTR, regulatory and splice-site classes need genomic input.
 *  An identifier whose stated reference residue differs from the transcript, or whose amino
 *  acid change needs more than one nucleotide change, is not annotated (vep warns and skips it).
 *
 *  Source of the transcripts (vep_mode): "database" queries the public Ensembl MySQL server
 *  (vep_db_host:vep_db_port), which needs internet but no download and is slow for many
 *  variants; "offline" reads a local cache; "auto" (the default) uses the database when the
 *  server answers and the offline cache otherwise.
 *
 *  Requires the `ensembl-vep` package in the task environment (the `vep` command, and
 *  `vep_install` for the offline cache). A process that does not find them writes the
 *  header-only table and says so in .command.err; install it with
 *    micromamba install -n phylophere -c conda-forge -c bioconda ensembl-vep
 *
 *  Offline cache: the VEP cache of a species and assembly is several GB (about 28 GB to
 *  download for human), so it lives in vep_cache_dir rather than in the work directory. An
 *  empty cache directory is populated once, on first use, with
 *    vep_install -a cf -s <species> -y <assembly> -c <vep_cache_dir> --NO_HTSLIB
 *  Point --vep_cache_dir at a populated cache to reuse or share one. When it is
 *  unset, the VEP workflow resolves it to ~/.cache/phylophere/vep/<species>_<assembly>.
 *
 *  Consumes:  position_scores.tsv (SCORING), directory of per-gene MAP files,
 *             gene_ensembl_file (gene, human_protein_id), VEP cache directory
 *  Produces:  vep/ensembl_vep_mapped.tsv (header only when ensembl-vep is not installed,
 *             the cache cannot be installed, vep fails or no HGVS identifier resolves)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Ensembl VEP annotation ───────────────────────────────────────────────────

process ENSEMBL_VEP_ANNOTATE {
    tag "ensembl_vep"
    label 'process_medium'
    errorStrategy 'ignore'   // a failed annotation never stops the run

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'ensembl_vep_mapped.tsv'

    input:
    path position_scores
    path vep_map_dir
    path gene_ensembl_file
    val vep_cache_dir

    output:
    path "ensembl_vep_mapped.tsv", emit: ensembl_vep_tsv

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    def species = params.vep_species ?: 'homo_sapiens'
    def assembly = params.vep_assembly ?: 'GRCh38'
    def vep_mode = params.vep_mode ?: 'auto'
    def db_host = params.vep_db_host ?: 'ensembldb.ensembl.org'
    def db_port = params.vep_db_port ?: 3306
    """
    cp ${local_dir}/build_vep_hgvs.py ${local_dir}/join_vep_output.py ${local_dir}/vep_common.py .

    header_only() {
        printf 'Gene\tPosition\tside\tUploaded_variation\tLocation\tAllele\tEnsembl_Gene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
    }

    # Without the ensembl-vep package there is nothing to run: say so before anything else.
    if ! command -v vep >/dev/null 2>&1; then
        echo "ERROR 'vep' is not on the PATH of this task: the ensembl-vep package is missing from the environment." >&2
        echo "      Install it with: micromamba install -n phylophere -c conda-forge -c bioconda ensembl-vep" >&2
        echo "      Skipping Ensembl VEP annotation; the output table has only its header." >&2
        header_only
        exit 0
    fi

    # Source of the transcripts: the public database when the server answers (auto, database), else the offline cache.
    mode="${vep_mode}"
    if [[ "\$mode" == "auto" ]]; then
        if timeout 15 bash -c "exec 3<>/dev/tcp/${db_host}/${db_port}" 2>/dev/null; then
            mode="database"
        else
            mode="offline"
        fi
        echo "INFO vep_mode auto resolved to \$mode (${db_host}:${db_port})." >&2
    fi

    if [[ "\$mode" == "database" ]]; then
        vep_source=(--database --host "${db_host}" --port ${db_port} --user anonymous)
    else
        if ! command -v vep_install >/dev/null 2>&1; then
            echo "ERROR 'vep_install' is not on the PATH of this task: the offline cache cannot be populated. Skipping Ensembl VEP annotation." >&2
            header_only
            exit 0
        fi
        mkdir -p "${vep_cache_dir}"
        if [[ -z "\$(ls -A "${vep_cache_dir}" 2>/dev/null)" ]]; then
            echo "INFO VEP cache empty at ${vep_cache_dir} -- populating once via vep_install (${species}/${assembly})." >&2
            if ! vep_install -a cf -s ${species} -y ${assembly} -c "${vep_cache_dir}" --NO_HTSLIB --NO_UPDATE --NO_TEST --QUIET; then
                echo "WARN vep_install failed for ${species}/${assembly} in ${vep_cache_dir} -- skipping Ensembl VEP annotation." >&2
                header_only
                exit 0
            fi
        fi
        vep_source=(--offline --cache --dir_cache "${vep_cache_dir}")
    fi

    # --canonical exists in every vep release; --mane from release 105. `vep --help` lists only the basic options, so
    # the release is read from its version block.
    vep_release="\$(vep --help 2>&1 | awk -F: '/ensembl-vep/ {gsub(/ /, "", \$2); print int(\$2); exit}')"
    mane_flag=()
    if [[ "\${vep_release:-0}" -ge 105 ]]; then mane_flag=(--mane); fi
    echo "INFO vep release \${vep_release:-unknown}, source \$mode." >&2

    python3 build_vep_hgvs.py "${position_scores}" "${vep_map_dir}" "${gene_ensembl_file}" hgvs_ids.txt id_map.tsv

    if [[ ! -s hgvs_ids.txt ]]; then
        echo "WARN No HGVS identifiers resolved (missing MAP files or protein IDs) — skipping Ensembl VEP call." >&2
        header_only
        exit 0
    fi

    if ! vep \
        --input_file hgvs_ids.txt \
        --format hgvs \
        --output_file vep_tab_output.txt \
        --tab \
        "\${vep_source[@]}" \
        --species ${species} \
        --assembly ${assembly} \
        --canonical \
        \${mane_flag[@]+"\${mane_flag[@]}"} \
        --force_overwrite \
        --no_stats; then
        echo "WARN vep failed -- the output table has only its header." >&2
        header_only
        exit 0
    fi

    python3 join_vep_output.py vep_tab_output.txt id_map.tsv ensembl_vep_mapped.tsv
    """

    stub:
    """
    printf 'Gene\tPosition\tside\tUploaded_variation\tLocation\tAllele\tEnsembl_Gene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
    """
}
