#!/usr/bin/env nextflow

/*
 * SCORING_REPORT
 * ──────────────
 * Render an HTML summary report for position-level and gene-level
 * CAAS scores.  Calls SCORING_report.Rmd.
 *
 * Inputs (15)
 * ──────
 *   position_scores      : path — position_scores.tsv
 *   gene_scores          : path — gene_scores.tsv
 *   gene_correlations    : path — gene_correlations.tsv
 *   stress_summary       : path — position_score_stress_summary.tsv (or NO_SCORING_STRESS_SUMMARY sentinel)
 *   stress_correlations  : path — position_score_stress_correlations.tsv (or NO_SCORING_STRESS_CORR sentinel)
 *   stress_rank_agreement: path — position_score_stress_rank_agreement.tsv (or NO_SCORING_STRESS_RANK sentinel)
 *   stress_top_overlap   : path — position_score_stress_top_overlap.tsv (or NO_SCORING_STRESS_OVERLAP sentinel)
 *   stress_variants      : path — position_score_stress_variants.tsv (or NO_SCORING_STRESS_VARIANTS sentinel)
 *   stress_latent_loadings: path — position_score_stress_latent_loadings.tsv (or NO_SCORING_STRESS_LOADINGS sentinel)
 *   fade_site_top_file   : path — per-site FADE BF TSV top direction (or NO_FADE_SITE_TOP sentinel)
 *   fade_site_bot_file   : path — per-site FADE BF TSV bottom direction (or NO_FADE_SITE_BOT sentinel)
 *   vep_transvar         : path — TransVar annotation TSV (or NO_VEP_TRANSVAR sentinel)
 *   vep_primateai        : path — PrimateAI-3D score TSV (or NO_VEP_PRIMATEAI sentinel)
 *   vep_aa2prot          : path — aa2prot_global.csv (or NO_VEP_AA2PROT sentinel)
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
    publishDir path: "${params.outdir}/HTML_reports",
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
    path stress_latent_loadings
    path fade_site_top_file  // optional: per-site FADE BF TSV top direction  (NO_FADE_SITE sentinel when absent)
    path fade_site_bot_file  // optional: per-site FADE BF TSV bottom direction (NO_FADE_SITE sentinel when absent)
    path vep_transvar        // optional: TransVar annotation TSV (NO_VEP_TRANSVAR sentinel when absent)
    path vep_primateai       // optional: PrimateAI-3D score TSV  (NO_VEP_PRIMATEAI sentinel when absent)
    path vep_aa2prot         // optional: aa2prot_global.csv — alignment→protein position map (NO_VEP_AA2PROT sentinel)
    path genomic_info        // optional: gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)

    output:
    path "SCORING_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def local_dir      = "${baseDir}/subworkflows/SCORING/local"
    def outdir         = "${params.outdir}/scoring"
    def traitname      = params.traitname ?: 'unknown_trait'
    def top_pct        = params.scoring_position_top_pct  ?: 0.10
    def top5_pct       = params.scoring_position_top5_pct ?: 0.05
    def top1_pct       = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct      = params.scoring_gene_top_pct      ?: 0.10
    def g_top5_pct     = params.scoring_gene_top5_pct     ?: 0.05
    def g_top1_pct     = params.scoring_gene_top1_pct     ?: 0.01
    // Resolve optional sentinel files: pass 'NULL' (R NULL) when no real file is staged
    def stress_summary_arg = (stress_summary.name =~ /^NO_SCORING_STRESS_SUMMARY/) ? 'NULL' : "'${stress_summary}'"
    def stress_corr_arg = (stress_correlations.name =~ /^NO_SCORING_STRESS_CORR/) ? 'NULL' : "'${stress_correlations}'"
    def stress_rank_arg = (stress_rank_agreement.name =~ /^NO_SCORING_STRESS_RANK/) ? 'NULL' : "'${stress_rank_agreement}'"
    def stress_overlap_arg = (stress_top_overlap.name =~ /^NO_SCORING_STRESS_OVERLAP/) ? 'NULL' : "'${stress_top_overlap}'"
    def stress_variants_arg = (stress_variants.name =~ /^NO_SCORING_STRESS_VARIANTS/) ? 'NULL' : "'${stress_variants}'"
    def stress_loadings_arg = (stress_latent_loadings.name =~ /^NO_SCORING_STRESS_LOADINGS/) ? 'NULL' : "'${stress_latent_loadings}'"
    def fs_top_arg = (fade_site_top_file.name =~ /^NO_FADE_SITE_TOP/) ? 'NULL' : "'${fade_site_top_file}'"
    def fs_bot_arg = (fade_site_bot_file.name =~ /^NO_FADE_SITE_BOT/) ? 'NULL' : "'${fade_site_bot_file}'"

    def tv_arg  = (vep_transvar.name  =~ /^NO_VEP_TRANSVAR/)  ? 'NULL' : "'${vep_transvar}'"
    def pai_arg = (vep_primateai.name =~ /^NO_VEP_PRIMATEAI/) ? 'NULL' : "'${vep_primateai}'"
    def a2p_arg = (vep_aa2prot.name   =~ /^NO_VEP_AA2PROT/)   ? 'NULL' : "'${vep_aa2prot}'"
    def gi_arg  = (genomic_info.name  =~ /^NO_GENOMIC_INFO/)  ? 'NULL' : "'${genomic_info}'"
    def win_size = params.scoring_window_size_bp ?: 1000000

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'SCORING_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    top5_pct             = ${top5_pct},
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top5_pct        = ${g_top5_pct},
                    gene_top1_pct        = ${g_top1_pct},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    stress_loadings_file = ${stress_loadings_arg},
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    vep_aa2prot_file     = ${a2p_arg},
                    genomic_info_file    = ${gi_arg},
                    window_size_bp       = ${win_size},
                    direction            = 'combined'
                ),
                output_file = 'SCORING_report_${traitname}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                'SCORING_report.Rmd',
                params = list(
                    position_scores_file = '${position_scores}',
                    gene_scores_file     = '${gene_scores}',
                    gene_corr_file       = '${gene_correlations}',
                    traitname            = '${traitname}',
                    output_dir           = '${outdir}',
                    top_pct              = ${top_pct},
                    top5_pct             = ${top5_pct},
                    top1_pct             = ${top1_pct},
                    gene_top_pct         = ${g_top_pct},
                    gene_top5_pct        = ${g_top5_pct},
                    gene_top1_pct        = ${g_top1_pct},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    stress_loadings_file = ${stress_loadings_arg},
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    vep_aa2prot_file     = ${a2p_arg},
                    genomic_info_file    = ${gi_arg},
                    window_size_bp       = ${win_size},
                    direction            = 'combined'
                ),
                output_file = 'SCORING_report_${traitname}.html'
            )
        "
        """
    }
}
