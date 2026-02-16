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
    def align_dir = params.ct_disambig_alignment_dir
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = params.ct_disambig_taxid_mapping ?: ''
    def ensembl_file = params.ct_disambig_ensembl_genes_file ?: ''
    def task_cpus = task.cpus ?: 1
    def threads = params.ct_disambig_threads ?: task_cpus
    def workers = params.ct_disambig_workers ?: task_cpus

    """
    mkdir -p ct_disambiguation
    cp -R ${local_dir}/* .

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
      ${asr_cache_dir ? "--asr-cache-dir ${asr_cache_dir}" : ''} \
      ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \
      ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}
