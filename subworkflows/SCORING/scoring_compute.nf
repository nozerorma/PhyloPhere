#!/usr/bin/env nextflow

/*
 * SCORING_COMPUTE
 * ───────────────
 * Compute composite CAAS scores at position-level and gene-level.
 * Integrates outputs from CT_POSTPROC, PGLS, FADE, RERConverge,
 * and CT_ACCUMULATION.
 *
 * Inputs
 * ──────
 *   postproc_file      : path — filtered_discovery.tsv (mandatory)
 *   pgls_file          : path — site_pgls.tsv (or NO_FILE sentinel)
 *   fade_summary_top   : path — fade_summary_top.tsv (or NO_FILE)
 *   fade_summary_bottom: path — fade_summary_bottom.tsv (or NO_FILE)
 *   rer_summary        : path — rerconverge_summary_{trait}.tsv (or NO_FILE)
 *   accum_files        : path — directory or collected CSVs (or NO_FILE)
 *
 * Outputs
 * ───────
 *   position_scores    : TSV with per Gene×Position scores
 *   gene_scores        : TSV with per Gene scores
 *   gene_correlations  : TSV with pairwise correlations
 *   position_gene_lists: TXT files for position-level ORA
 *   gene_gene_lists    : TXT files for gene-level ORA
 */

process SCORING_COMPUTE {
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}"
    label 'error_retry'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/scoring/gene_lists/position",
               mode: 'copy', overwrite: true,
               pattern: 'pos_gene_lists/*.txt'
    publishDir path: "${params.outdir}/scoring/gene_lists/gene",
               mode: 'copy', overwrite: true,
               pattern: 'gene_gene_lists/*.txt'

    input:
    path postproc_file
    path pgls_file
    path fade_summary_top
    path fade_summary_bottom
    path rer_summary
    path accum_files

    output:
    path "position_scores.tsv",    emit: position_scores
    path "gene_scores.tsv",        emit: gene_scores
    path "gene_correlations.tsv",  emit: gene_correlations
    path "pos_gene_lists/*.txt",   emit: position_gene_lists, optional: true
    path "gene_gene_lists/*.txt",  emit: gene_gene_lists,     optional: true

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local"
    def top_pct    = params.scoring_position_top_pct  ?: 0.10
    def top1_pct   = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct  = params.scoring_gene_top_pct      ?: 0.10
    def g_top1_pct = params.scoring_gene_top1_pct     ?: 0.01
    // Determine accum_dir: if accum_files is a real directory, pass it;
    // otherwise check if multiple CSVs were staged (Nextflow stages them flat).
    // accum_files may be a single sentinel file or a list of CSVs
    def accum_arg  = (accum_files instanceof List
                        ? (accum_files.size() == 1 && accum_files[0].name.startsWith('NO_') ? 'NO_ACCUM' : '.')
                        : (accum_files.name.startsWith('NO_') ? 'NO_ACCUM' : '.'))

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/scoring_compute.R .

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc   '${postproc_file}' \
            --pgls       '${pgls_file}' \
            --fade_top   '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer        '${rer_summary}' \
            --accum_dir  '${accum_arg}' \
            --top_pct    ${top_pct} \
            --top1_pct   ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top1_pct ${g_top1_pct}
        """
    } else {
        """
        cp ${local_dir}/scoring_compute.R .

        Rscript scoring_compute.R \
            --postproc   '${postproc_file}' \
            --pgls       '${pgls_file}' \
            --fade_top   '${fade_summary_top}' \
            --fade_bottom '${fade_summary_bottom}' \
            --rer        '${rer_summary}' \
            --accum_dir  '${accum_arg}' \
            --top_pct    ${top_pct} \
            --top1_pct   ${top1_pct} \
            --gene_top_pct  ${g_top_pct} \
            --gene_top1_pct ${g_top1_pct}
        """
    }
}
