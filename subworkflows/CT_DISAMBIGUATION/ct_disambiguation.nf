#!/usr/bin/env nextflow

/*
 * CT disambiguation subworkflow
 */

process CT_DISAMBIGUATION_RUN {
    tag "ct_disambiguation"
    label 'process_resample'
    publishDir path: "${params.outdir}/ct_disambiguation", mode: 'copy', overwrite: true

    input:
    path meta_caas
    path trait_file
    path tree_file

    output:
    path("ct_disambiguation"), emit: results_dir
    path("ct_disambiguation/caas_convergence_master.csv"), emit: master_csv

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def disambig_script = params.ct_disambig_script ?: "${local_dir}/disambiguation_main.py"
    def align_dir = params.alignment
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''

    def asr_mode = params.ct_disambig_asr_mode
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def task_cpus = task.cpus ?: 1
    def threads = params.ct_disambig_threads ?: task_cpus
    def workers = params.ct_disambig_workers ?: task_cpus

    """
    mkdir -p ct_disambiguation
    cp -R ${local_dir}/* .
    # Remove stale .pyc / __pycache__ dirs so Python always compiles from source
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    echo "[ct_disambiguation] Inputs:"
    echo "  meta_caas=${meta_caas}"
    echo "  trait_file=${trait_file}"
    echo "  tree_file=${tree_file}"

    # Validate ASR mode / cache dir combination
    if [ -z "${asr_cache_dir}" ]; then
      echo "ERROR: ct_disambig_asr_cache_dir must be set (current asr_mode: '${params.ct_disambig_asr_mode}')" >&2
      exit 1
    fi

    if [ ! -s "${meta_caas}" ]; then
      echo "ERROR: metadata file is missing or empty: ${meta_caas}" >&2
      exit 1
    fi

    if [ ! -s "${trait_file}" ]; then
      echo "ERROR: trait file is missing or empty: ${trait_file}" >&2
      exit 1
    fi

    meta_rows=\$(wc -l < "${meta_caas}")
    trait_rows=\$(wc -l < "${trait_file}")

    echo "  meta_rows=\${meta_rows}"
    echo "  trait_rows=\${trait_rows}"

    if [ "\${meta_rows}" -le 1 ]; then
      echo "ERROR: metadata file has header only (no data rows): ${meta_caas}" >&2
      exit 1
    fi

    if [ "\${trait_rows}" -le 1 ]; then
      echo "ERROR: trait file has header only (no data rows): ${trait_file}" >&2
      exit 1
    fi

    python3 ./disambiguation_main.py \
      --alignment-dir ${align_dir} \
      --tree ${tree_file} \
      --caas-metadata ${meta_caas} \
      --trait-file ${trait_file} \
      --output-dir ct_disambiguation \
      --asr-mode ${params.ct_disambig_asr_mode} \
      --asr-model ${params.ct_disambig_asr_model} \
      --convergence-mode ${params.ct_disambig_convergence_mode} \
      --posterior-threshold ${params.ct_disambig_posterior_threshold} \
      --threads ${threads} \
      --workers ${workers} \
      --max-tasks-per-child ${params.ct_disambig_max_tasks_per_child} \
      ${params.ct_disambig_include_non_significant ? '--include-non-significant' : ''} \
      ${params.ct_disambig_run_diagnostics ? '--run-diagnostics' : ''} \
      ${params.ct_disambig_skip_gene_lists ? '--skip-gene-lists' : ''} \
      ${params.ct_disambig_verbose ? '--verbose' : ''} \
      ${params.ct_disambig_use_all_mrca_filter ? '--use-all-mrca-filter' : ''} \
      ${asr_cache_dir ? "--asr-cache-dir ${asr_cache_dir}" : ''} \
      ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \
      ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}
