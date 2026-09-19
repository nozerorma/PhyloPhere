#!/usr/bin/env nextflow

/*
 * ENSEMBL_VEP_ANNOTATE
 * ────────────────────
 * Annotates CAAS ancestral->derived amino-acid changes with Ensembl VEP's own
 * consequence prediction — independent of PrimateAI-3D/COSMIC, so consequence
 * annotation is available even when neither pathogenicity database is
 * supplied. Not to be confused with this pipeline's own --vep toggle (which
 * gates this whole workflow); "Ensembl VEP" here means the official
 * variant_effect_predictor CLI tool.
 *
 * Unlike PRIMATEAI_MAP/COSMIC_MAP, the ancestral (reference) amino acid comes
 * from this pipeline's own ASR descriptor (build_vep_hgvs.py), not from an
 * external pathogenicity database's bundled reference annotation — so this
 * doesn't introduce a new external reference-proteome dependency.
 *
 * Requires a local, pre-downloaded VEP cache (--vep_cache_dir): VEP's offline
 * per-species/assembly cache is multi-GB, the same "cache large, don't commit"
 * pattern used for STRING/eggNOG/Pfam-A elsewhere in this pipeline. Populate
 * it once with:
 *   vep_install -a cf -s <species> -y <assembly> -c <vep_cache_dir> --NO_HTSLIB
 */

process ENSEMBL_VEP_ANNOTATE {
    tag "ensembl_vep"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/vep",
               mode: 'copy', overwrite: true,
               pattern: 'ensembl_vep_mapped.tsv'

    input:
    path caas_file
    path vep_map_dir
    path gene_ensembl_file
    path vep_cache_dir

    output:
    path "ensembl_vep_mapped.tsv", emit: ensembl_vep_tsv

    stub:
    """
    printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/VEP/local/src"
    def species = params.vep_species ?: 'homo_sapiens'
    def assembly = params.vep_assembly ?: 'GRCh38'
    """
    cp ${local_dir}/build_vep_hgvs.py ${local_dir}/join_vep_output.py ${local_dir}/vep_common.py .

    if [[ ! -d "${vep_cache_dir}" || -z "\$(ls -A "${vep_cache_dir}" 2>/dev/null)" ]]; then
        echo "WARN Missing/empty VEP cache directory: ${vep_cache_dir}. Skipping Ensembl VEP annotation." >&2
        echo "WARN Populate it once with: vep_install -a cf -s ${species} -y ${assembly} -c ${vep_cache_dir} --NO_HTSLIB" >&2
        printf 'Gene\tPosition\tcaap_group\tUploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n' > ensembl_vep_mapped.tsv
        exit 0
    fi

    python3 build_vep_hgvs.py "${caas_file}" "${vep_map_dir}" "${gene_ensembl_file}" hgvs_ids.txt id_map.tsv

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
}
