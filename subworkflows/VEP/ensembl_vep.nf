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
 *  Steps: build_vep_hgvs.py writes one protein-level HGVS identifier per change
 *  (the ancestral and derived residues are the ASR ones of position_scores.tsv, so no
 *  reference proteome is needed), `vep --offline` annotates them, and join_vep_output.py
 *  attaches the result to Gene, Position and caap_group.
 *
 *  Requires the `ensembl-vep` package (the `vep` and `vep_install` commands) in the task
 *  environment. A process that does not find them writes the header-only table and says
 *  so in .command.err; install it with
 *    micromamba install -n phylophere -c conda-forge -c bioconda ensembl-vep
 *
 *  Offline cache: the VEP cache of a species and assembly is several GB, so it
 *  lives in vep_cache_dir rather than in the work directory. An empty cache
 *  directory is populated once, on first use, with
 *    vep_install -a cf -s <species> -y <assembly> -c <vep_cache_dir> --NO_HTSLIB
 *  Point --vep_cache_dir at a populated cache to reuse or share one. When it is
 *  unset, the VEP workflow resolves it to ~/.cache/phylophere/vep/<species>_<assembly>.
 *
 *  Consumes:  position_scores.tsv (SCORING), directory of per-gene MAP files,
 *             gene_ensembl_file (gene, human_protein_id), VEP cache directory
 *  Produces:  vep/ensembl_vep_mapped.tsv (header only when ensembl-vep is not installed,
 *             the cache cannot be installed or no HGVS identifier resolves)
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
    """
    cp ${local_dir}/build_vep_hgvs.py ${local_dir}/join_vep_output.py ${local_dir}/vep_common.py .

    # Without the ensembl-vep package there is nothing to run: say so before anything else.
    for tool in vep vep_install; do
        if ! command -v "\$tool" >/dev/null 2>&1; then
            echo "ERROR '\$tool' is not on the PATH of this task: the ensembl-vep package is missing from the environment." >&2
            echo "      Install it with: micromamba install -n phylophere -c conda-forge -c bioconda ensembl-vep" >&2
            echo "      Skipping Ensembl VEP annotation; the output table has only its header." >&2
            printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
            exit 0
        fi
    done

    mkdir -p "${vep_cache_dir}"
    if [[ -z "\$(ls -A "${vep_cache_dir}" 2>/dev/null)" ]]; then
        echo "INFO VEP cache empty at ${vep_cache_dir} -- populating once via vep_install (${species}/${assembly})." >&2
        if ! vep_install -a cf -s ${species} -y ${assembly} -c "${vep_cache_dir}" --NO_HTSLIB --NO_UPDATE --NO_TEST --QUIET; then
            echo "WARN vep_install failed for ${species}/${assembly} in ${vep_cache_dir} -- skipping Ensembl VEP annotation." >&2
            printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
            exit 0
        fi
    fi

    python3 build_vep_hgvs.py "${position_scores}" "${vep_map_dir}" "${gene_ensembl_file}" hgvs_ids.txt id_map.tsv

    if [[ ! -s hgvs_ids.txt ]]; then
        echo "WARN No HGVS identifiers resolved (missing MAP files or protein IDs) — skipping Ensembl VEP call." >&2
        printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
        exit 0
    fi

    vep \
        --input_file hgvs_ids.txt \
        --format hgvs \
        --output_file vep_tab_output.txt \
        --tab \
        --offline \
        --cache \
        --dir_cache "${vep_cache_dir}" \
        --species ${species} \
        --assembly ${assembly} \
        --force_overwrite \
        --no_stats

    python3 join_vep_output.py vep_tab_output.txt id_map.tsv ensembl_vep_mapped.tsv
    """

    stub:
    """
    printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
    """
}
