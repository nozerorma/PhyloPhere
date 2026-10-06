#!/usr/bin/env nextflow

// ct_evidence.nf — Evidence of the best-scored positions: what each domain of each hypothesis saw.
// PhyloPhere | subworkflows/CT_DISAMBIGUATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CAAS_EVIDENCE: for the N best positions of position_scores.tsv, re-scores their rows
 *  of the observed discovery.tab with the code of the observed labeling
 *  (explain_positions.py), keeping the rows before the hypotheses of a position are
 *  pooled.
 *
 *  Runs after SCORING, only when params.caas_evidence_top_n is greater than 0 (main.nf).
 *
 *  Consumes:  discovery.tab, scoring/position_scores.tsv, observed design (trait file or
 *             traitfile_H*.tab directory), species tree
 *  Produces:  evidence/ with evidence_top<N>.tsv (one row per entry and domain) and
 *             top_positions.tsv, published to <outdir>/scoring/evidence
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Evidence of the top positions ────────────────────────────────────────────

process CAAS_EVIDENCE {
    tag "caas_evidence"
    label 'process_low'
    // The directory is published by name: a `dir/**` pattern publishes nothing for a process that declares only the directory.
    publishDir path: "${params.outdir}/scoring", mode: 'copy', overwrite: true, pattern: 'evidence'

    input:
    path discovery          // the observed discovery.tab
    path position_scores    // scoring/position_scores.tsv
    path design             // observed design: trait file, or the directory of traitfile_H*.tab
    path tree_file

    output:
    path "evidence",                  emit: evidence_dir
    path "evidence/evidence_top*.tsv", emit: evidence_tsv
    path "evidence/top_positions.tsv", emit: top_positions

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    # One worker per gene; BLAS and OpenMP stay at one thread each.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cp -R ${local_dir}/* .
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    if [ -z "${asr_cache_dir}" ]; then
      echo "ERROR: ct_disambig_asr_cache_dir must be set" >&2
      exit 1
    fi

    ${run} ./explain_positions.py \\
        --alignment-dir ${params.alignment} \\
        --tree ${tree_file} \\
        --discovery ${discovery} \\
        --position-scores ${position_scores} \\
        --top ${params.caas_evidence_top_n as int} \\
        --design ${design} \\
        --output-dir evidence \\
        --asr-model ${params.ct_disambig_asr_model} \\
        --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
        --workers ${task.cpus} \\
        --asr-cache-dir ${asr_cache_dir} \\
        ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
        ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    """
}
