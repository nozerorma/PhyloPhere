#!/usr/bin/env nextflow

/*
 * SCORING_REPORT
 * ──────────────
 * Render an HTML summary report for position-level and gene-level
 * CAAS scores.  Calls 11.Scoring_report.Rmd.
 *
 * Inputs
 * ──────
 *   position_scores      : path — position_scores.tsv
 *   gene_scores          : path — gene_scores.tsv
 *   gene_correlations    : path — gene_correlations.tsv
 *   stress_summary       : path — position_score_stress_summary.tsv (or NO_SCORING_STRESS_SUMMARY sentinel)
 *   stress_correlations  : path — position_score_stress_correlations.tsv (or NO_SCORING_STRESS_CORR sentinel)
 *   stress_rank_agreement: path — position_score_stress_rank_agreement.tsv (or NO_SCORING_STRESS_RANK sentinel)
 *   stress_top_overlap   : path — position_score_stress_top_overlap.tsv (or NO_SCORING_STRESS_OVERLAP sentinel)
 *   stress_variants      : path — position_score_stress_variants.tsv (or NO_SCORING_STRESS_VARIANTS sentinel)
 *   fade_site_top_file   : path — per-site FADE BF TSV top direction (or NO_FADE_SITE_TOP sentinel)
 *   fade_site_bot_file   : path — per-site FADE BF TSV bottom direction (or NO_FADE_SITE_BOT sentinel)
 *   genomic_info         : path — gene genomic coords TSV (or NO_GENOMIC_INFO sentinel)
 *
 * Outputs
 * ───────
 *   report : HTML report
 */

process SCORING_REPORT {
    tag "scoring_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path position_scores
    path gene_scores
    path gene_correlations
    path stress_summary
    path stress_correlations
    path stress_rank_agreement
    path stress_top_overlap
    path stress_variants
    path fade_site_top_file  // optional: per-site FADE BF TSV top direction  (NO_FADE_SITE sentinel when absent)
    path fade_site_bot_file  // optional: per-site FADE BF TSV bottom direction (NO_FADE_SITE sentinel when absent)
    path genomic_info        // optional: gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)
    path caas_perms          // optional: CAAS permulation RDS (asr + caas null) (NO_FILE/NO_CAAS_PERMS sentinel when absent)
    path caas_pos_pval       // optional: LOO null_pvalue_boot per (gene,position,scheme) (NO_CAAS_POS_PVAL sentinel)
    path caas_pos_sample     // optional: cycle-stratified per-scheme sample for distribution plots (NO_CAAS_POS_SAMPLE sentinel)
    path caas_pos_quantiles  // optional: per (cycle,scheme) null distribution shape (NO_CAAS_POS_QUANTILES sentinel)
    path filtered_discovery  // observed per-(gene,position,scheme) asr_path_score (filtered_discovery.tsv) for the null overlay
    path background_file

    output:
    path "11.Scoring_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def local_dir      = "${baseDir}/subworkflows/SCORING/local"
    def outdir         = "${params.outdir}/scoring"
    def traitname      = params.traitname ?: 'unknown_trait'
    def top_pct        = params.scoring_position_top_pct   ?: 0.10
    def g_top_pct      = params.scoring_gene_top_pct       ?: 0.10
    // Resolve optional sentinel files: pass 'NULL' (R NULL) when no real file is staged
    def stress_summary_arg = (stress_summary.name =~ /^NO_SCORING_STRESS_SUMMARY/) ? 'NULL' : "'${stress_summary}'"
    def stress_corr_arg = (stress_correlations.name =~ /^NO_SCORING_STRESS_CORR/) ? 'NULL' : "'${stress_correlations}'"
    def stress_rank_arg = (stress_rank_agreement.name =~ /^NO_SCORING_STRESS_RANK/) ? 'NULL' : "'${stress_rank_agreement}'"
    def stress_overlap_arg = (stress_top_overlap.name =~ /^NO_SCORING_STRESS_OVERLAP/) ? 'NULL' : "'${stress_top_overlap}'"
    def stress_variants_arg = (stress_variants.name =~ /^NO_SCORING_STRESS_VARIANTS/) ? 'NULL' : "'${stress_variants}'"
    def fs_top_arg = (fade_site_top_file.name =~ /^NO_FADE_SITE_TOP/) ? 'NULL' : "'${fade_site_top_file}'"
    def fs_bot_arg = (fade_site_bot_file.name =~ /^NO_FADE_SITE_BOT/) ? 'NULL' : "'${fade_site_bot_file}'"

    def gi_arg  = (genomic_info.name  =~ /^NO_GENOMIC_INFO/)  ? 'NULL' : "'${genomic_info}'"
    def perms_arg = (caas_perms.name =~ /^NO_CAAS_PERMS|^NO_FILE/) ? 'NULL' : "'${caas_perms}'"
    def pos_pval_arg = (caas_pos_pval.name =~ /^NO_CAAS_POS_PVAL|^NO_FILE/) ? 'NULL' : "'${caas_pos_pval}'"
    def pos_sample_arg = (caas_pos_sample.name =~ /^NO_CAAS_POS_SAMPLE|^NO_FILE/) ? 'NULL' : "'${caas_pos_sample}'"
    def pos_quantiles_arg = (caas_pos_quantiles.name =~ /^NO_CAAS_POS_QUANTILES|^NO_FILE/) ? 'NULL' : "'${caas_pos_quantiles}'"
    def filt_disc_arg = (filtered_discovery.name =~ /^NO_POSTPROC|^NO_FILE/) ? 'NULL' : "'${filtered_discovery}'"
    def bg_file_arg = (background_file.name =~ /^NO_BACKGROUND|^NO_FILE/) ? 'NULL' : "'${background_file}'"
    def win_size = params.scoring_window_size_bp ?: 1000000

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '11.Scoring_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    gene_top_pct         = ${g_top_pct},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    genomic_info_file    = ${gi_arg},
                    caas_perms_file      = ${perms_arg},
                    caas_pos_pval_file   = ${pos_pval_arg},
                    caas_pos_sample_file = ${pos_sample_arg},
                    caas_pos_quantiles_file = ${pos_quantiles_arg},
                    filtered_discovery_file = ${filt_disc_arg},
                    background_file      = ${bg_file_arg},
                    window_size_bp       = ${win_size},
                    direction            = 'combined',
                    seed                 = '${params.seed ?: 1998}'
                ),
                output_file = '11.Scoring_report_${traitname}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                '11.Scoring_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    gene_top_pct         = ${g_top_pct},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    genomic_info_file    = ${gi_arg},
                    caas_perms_file      = ${perms_arg},
                    caas_pos_pval_file   = ${pos_pval_arg},
                    caas_pos_sample_file = ${pos_sample_arg},
                    caas_pos_quantiles_file = ${pos_quantiles_arg},
                    filtered_discovery_file = ${filt_disc_arg},
                    background_file      = ${bg_file_arg},
                    window_size_bp       = ${win_size},
                    direction            = 'combined',
                    seed                 = '${params.seed ?: 1998}'
                ),
                output_file = '11.Scoring_report_${traitname}.html'
            )
        "
        """
    }
}
