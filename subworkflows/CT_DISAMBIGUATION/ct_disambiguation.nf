#!/usr/bin/env nextflow

/*
 * CT disambiguation subworkflow
 */

process CT_DISAMBIGUATION_RUN {
    tag "ct_disambiguation"
    label 'process_resample'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path meta_caas
    path trait_file
    path tree_file
    path hyp_pairs

    output:
    path("ct_disambiguation"), emit: results_dir
    path("ct_disambiguation/caas_convergence_master.csv"), emit: master_csv

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def disambig_script = "${local_dir}/disambiguation_main.py"
    def align_dir = params.alignment
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''

    def asr_mode = params.ct_disambig_asr_mode
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def task_cpus = task.cpus ?: 1
    def threads = task_cpus
    def workers = task_cpus

    """
    # Disambiguation is an irregular pure-Python tree walk over dict posteriors
    # (no Level-3 BLAS) and fans genes across an mp.Pool of ${workers} workers.
    # Pin math-library threads to 1 so a per-worker BLAS/OpenMP bloom can't
    # multiply by the worker count. Per-stage only (RERConverge needs multithread
    # BLAS elsewhere and must not see a global pin).
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
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

    if [ ! -e "${trait_file}" ]; then
      echo "ERROR: trait file or directory is missing: ${trait_file}" >&2
      exit 1
    fi

    meta_rows=\$(wc -l < "${meta_caas}")
    if [ -d "${trait_file}" ]; then
      trait_rows=\$(find -L "${trait_file}" -maxdepth 1 -type f -name '*.tab' -exec wc -l {} + | awk 'END {print \$1+0}')
      echo "  meta_rows=\${meta_rows}"
      echo "  trait_rows=\${trait_rows} (across all traitfiles in directory ${trait_file})"
    else
      trait_rows=\$(wc -l < "${trait_file}")
      echo "  meta_rows=\${meta_rows}"
      echo "  trait_rows=\${trait_rows}"
    fi

    if [ "\${meta_rows}" -le 1 ]; then
      echo "ERROR: metadata file has header only (no data rows): ${meta_caas}" >&2
      exit 1
    fi

    if [ "\${trait_rows}" -le 0 ]; then
      echo "ERROR: trait file/directory has no data rows: ${trait_file}" >&2
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
      --posterior-threshold ${params.ct_disambig_posterior_threshold} \
      --threads ${threads} \
      --workers ${workers} \
      --max-tasks-per-child ${params.ct_disambig_max_tasks_per_child} \
      --run-diagnostics \
      --verbose \
      ${hyp_pairs.name.startsWith('NO_') ? '' : "--hypotheses-pairs ${hyp_pairs}"} \
      ${asr_cache_dir ? "--asr-cache-dir ${asr_cache_dir}" : ''} \
      ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \
      ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}

// ── Batched by gene: split the metadata table into N gene subsets and run one
//    CT_DISAMBIGUATION_RUN-equivalent Nextflow task per subset, so one bad
//    batch (OOM, a pathological gene) retries on its own instead of
//    re-running the whole ~16k-gene sweep. See CT_DISAMBIGUATION_SPLIT_GENES
//    for how the metadata table is partitioned and CT_DISAMBIGUATION_MERGE
//    for how the batches' outputs are recombined.
process CT_DISAMBIGUATION_SPLIT_GENES {
    tag "split_by_gene (batch_size=${batchSize})"
    label 'process_low'

    input:
    path meta_caas
    val batchSize

    output:
    path("batch_*.meta_caas.tsv"), emit: batches

    script:
    def scripts_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/scripts"
    """
    python3 ${scripts_dir}/split_meta_caas_by_genes.py \
      --meta-caas ${meta_caas} \
      --batch-size ${batchSize} \
      --outdir .
    """
}

process CT_DISAMBIGUATION_RUN_BATCHED {
    tag "$meta_caas_batch"
    label 'process_resample'

    input:
    path meta_caas_batch
    path trait_file
    path tree_file
    path hyp_pairs

    output:
    path("ct_disambiguation"), emit: results_dir

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def align_dir = params.alignment
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''

    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def task_cpus = task.cpus ?: 1
    def threads = task_cpus
    def workers = task_cpus

    """
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    mkdir -p ct_disambiguation
    cp -R ${local_dir}/* .
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    if [ -z "${asr_cache_dir}" ]; then
      echo "ERROR: ct_disambig_asr_cache_dir must be set (current asr_mode: '${params.ct_disambig_asr_mode}')" >&2
      exit 1
    fi

    python3 ./disambiguation_main.py \
      --alignment-dir ${align_dir} \
      --tree ${tree_file} \
      --caas-metadata ${meta_caas_batch} \
      --trait-file ${trait_file} \
      --output-dir ct_disambiguation \
      --asr-mode ${params.ct_disambig_asr_mode} \
      --asr-model ${params.ct_disambig_asr_model} \
      --posterior-threshold ${params.ct_disambig_posterior_threshold} \
      --threads ${threads} \
      --workers ${workers} \
      --max-tasks-per-child ${params.ct_disambig_max_tasks_per_child} \
      --run-diagnostics \
      --verbose \
      ${hyp_pairs.name.startsWith('NO_') ? '' : "--hypotheses-pairs ${hyp_pairs}"} \
      ${asr_cache_dir ? "--asr-cache-dir ${asr_cache_dir}" : ''} \
      ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \
      ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}

// Batches partition genes disjointly (CT_DISAMBIGUATION_SPLIT_GENES assigns
// each gene to exactly one batch), so every piece here is a plain row-concat
// or directory union -- see merge_disambiguation_batches.py's module
// docstring for why no cross-gene aggregation is needed (unlike
// CAAS_PERMS_MERGE_DETAIL's downstream CAAS_PERMS_REBUILD step).
process CT_DISAMBIGUATION_MERGE {
    tag "ct_disambiguation_merge"
    label 'process_medium'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path batchDirs, stageAs: 'batch_*'

    output:
    path("ct_disambiguation"), emit: results_dir
    path("ct_disambiguation/caas_convergence_master.csv"), emit: master_csv

    script:
    def scripts_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/scripts"
    """
    python3 ${scripts_dir}/merge_disambiguation_batches.py \
      --batch-dirs batch_*/ \
      --output-dir ct_disambiguation
    """
}

// Regenerates plots/ once against the merged master CSV. Each batch's own
// disambiguation_main.py run already generated a (batch-local, incomplete)
// plots/ directory as a side effect of CT_DISAMBIGUATION_RUN_BATCHED --
// CT_DISAMBIGUATION_MERGE does not carry those through, so this replaces
// them with one run over the full, merged data.
process CT_DISAMBIGUATION_PLOTS {
    tag "ct_disambiguation_plots"
    label 'process_reporting'
    publishDir path: "${params.outdir}/ct_disambiguation", mode: 'copy', overwrite: true, pattern: 'plots/**'

    input:
    path merged_dir

    output:
    path("plots"), emit: plots, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scripts_dir = "${local_dir}/scripts"
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    """
    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    python3 ${scripts_dir}/regenerate_disambiguation_plots.py \
      --caas-csv ${merged_dir}/caas_convergence_master.csv \
      --output-dir . \
      ${asr_cache_dir ? "--asr-cache-dir ${asr_cache_dir}" : ''} \
      --node-dumps-root ${merged_dir}/diagnostics/node_dumps \
      ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}
