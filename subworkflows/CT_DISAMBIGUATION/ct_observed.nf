#!/usr/bin/env nextflow
// ct_observed.nf — Score the observed labeling of a discovery.tab that already exists.
// PhyloPhere | subworkflows/CT_DISAMBIGUATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CAAS_OBSERVED: scores the observed labeling of a discovery.tab given by the user.
 *
 *  A live run gets its observed results as the b_0 slice of the permulation core
 *  (CAAS_CORE_BATCHED and CAAS_CORE_MERGE). A run that reuses a discovery.tab has no
 *  replay, so its rows are scored here by the same code (observed_b0_main.py) and the
 *  master and meta_caas tables are written by the same adapters (contract_main.py).
 *
 *  Consumes:  discovery.tab, observed design (trait file or traitfile_H*.tab directory),
 *             species tree, hypothesis pairs (or a NO_* sentinel)
 *  Produces:  ct_disambiguation/ (master CSV), meta_caas/ (global_meta_caas.tsv)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Observed scoring ───────────────────────────────────────────────────────────

process CAAS_OBSERVED {
    tag "caas_observed"
    label 'process_resample'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, pattern: 'ct_disambiguation'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: 'meta_caas'

    input:
    path discovery
    path trait_file    // observed design: trait file, or the directory of traitfile_H*.tab
    path tree_file
    path hyp_pairs     // contrast_hypotheses_pairs.tsv or a NO_* sentinel (equal-weight node pooling)
    path taxid_map     // tax_id map of the species (curated by NAME_CURATION, or params.tax_id) or the NO_FILE sentinel

    output:
    path "ct_disambiguation", emit: results_dir
    path "ct_disambiguation/caas_convergence_master.csv", emit: master_csv
    path "meta_caas", emit: meta_caas
    path "meta_caas/global_meta_caas.tsv", emit: global_meta_caas

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def taxid_mapping = taxid_map.name != 'NO_FILE' ? taxid_map : ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    # One worker per gene: BLAS and OpenMP stay at one thread each so the workers do not oversubscribe the cpus.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cp -R ${local_dir}/* .
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    if [ -z "${asr_cache_dir}" ]; then
      echo "ERROR: ct_disambig_asr_cache_dir must be set" >&2
      exit 1
    fi
    mkdir -p "${asr_cache_dir}"
    if [ "\$(wc -l < "${discovery}")" -le 1 ]; then
      echo "ERROR: the discovery file has no data rows: ${discovery}" >&2
      exit 1
    fi
    if [ ! -e "${trait_file}" ]; then
      echo "ERROR: trait file or directory is missing: ${trait_file}" >&2
      exit 1
    fi

    ${run} ./observed_b0_main.py \\
        --alignment-dir ${params.alignment} \\
        --tree ${tree_file} \\
        --discovery ${discovery} \\
        --design ${trait_file} \\
        --output-dir shards \\
        --asr-model ${params.ct_disambig_asr_model} \\
        --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
        --workers ${task.cpus} \\
        --max-tasks-per-child ${params.ct_disambig_max_tasks_per_child} \\
        --asr-cache-dir ${asr_cache_dir} \\
        ${hyp_pairs.name.startsWith('NO_') ? '' : "--fop-pairs ${hyp_pairs}"} \\
        ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
        ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}

    ${run} ./contract_main.py --b0-dirs shards --design ${trait_file} --discovery-file ${discovery} --output-dir .
    """
}
