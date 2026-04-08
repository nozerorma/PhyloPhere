#!/usr/bin/env nextflow

/*
 * SCORING_REPORT
 * ──────────────
 * Render an HTML summary report for position-level and gene-level
 * CAAS scores.  Calls SCORING_report.Rmd.
 *
 * Inputs
 * ──────
 *   position_scores  : path — position_scores.tsv
 *   gene_scores      : path — gene_scores.tsv
 *   gene_correlations: path — gene_correlations.tsv
 *
 * Outputs
 * ───────
 *   report : HTML report
 */

process SCORING_REPORT {
    tag "scoring_report|${params.traitname ?: 'unknown_trait'}|${direction}"
    label 'process_reporting'

    publishDir path: { "${params.outdir}/scoring/${direction}" },
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    val  direction
    path position_scores
    path gene_scores
    path gene_correlations
    path pgls_excess_report  // optional: pgls_excess_for_report.tsv (NO_PGLS_EXCESS sentinel when absent)
    path stress_summary
    path stress_correlations
    path stress_rank_agreement
    path stress_top_overlap
    path stress_variants
    path stress_latent_loadings
    path vep_transvar        // optional: TransVar annotation TSV (NO_VEP_TRANSVAR sentinel when absent)
    path vep_primateai       // optional: PrimateAI-3D score TSV  (NO_VEP_PRIMATEAI sentinel when absent)
    path vep_aa2prot         // optional: aa2prot_global.csv — alignment→protein position map (NO_VEP_AA2PROT sentinel)
    path genomic_info        // optional: gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)

    output:
    path "SCORING_report_${direction}.html", emit: report

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
    def exc_arg = (pgls_excess_report.name =~ /^NO_PGLS_EXCESS/) ? 'NULL' : "'${pgls_excess_report}'"
    def stress_summary_arg = (stress_summary.name =~ /^NO_SCORING_STRESS_SUMMARY/) ? 'NULL' : "'${stress_summary}'"
    def stress_corr_arg = (stress_correlations.name =~ /^NO_SCORING_STRESS_CORR/) ? 'NULL' : "'${stress_correlations}'"
    def stress_rank_arg = (stress_rank_agreement.name =~ /^NO_SCORING_STRESS_RANK/) ? 'NULL' : "'${stress_rank_agreement}'"
    def stress_overlap_arg = (stress_top_overlap.name =~ /^NO_SCORING_STRESS_OVERLAP/) ? 'NULL' : "'${stress_top_overlap}'"
    def stress_variants_arg = (stress_variants.name =~ /^NO_SCORING_STRESS_VARIANTS/) ? 'NULL' : "'${stress_variants}'"
    def stress_loadings_arg = (stress_latent_loadings.name =~ /^NO_SCORING_STRESS_LOADINGS/) ? 'NULL' : "'${stress_latent_loadings}'"
    def tv_arg  = (vep_transvar.name  =~ /^NO_VEP_TRANSVAR/)  ? 'NULL' : "'${vep_transvar}'"
    def pai_arg = (vep_primateai.name =~ /^NO_VEP_PRIMATEAI/) ? 'NULL' : "'${vep_primateai}'"
    def a2p_arg = (vep_aa2prot.name   =~ /^NO_VEP_AA2PROT/)   ? 'NULL' : "'${vep_aa2prot}'"
    def gi_arg  = (genomic_info.name  =~ /^NO_GENOMIC_INFO/)  ? 'NULL' : "'${genomic_info}'"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
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
                    pgls_excess_file     = ${exc_arg},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    stress_loadings_file = ${stress_loadings_arg},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    vep_aa2prot_file     = ${a2p_arg},
                    genomic_info_file    = ${gi_arg},
                    direction            = '${direction}'
                ),
                output_file = 'SCORING_report_${direction}.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
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
                    pgls_excess_file     = ${exc_arg},
                    stress_summary_file  = ${stress_summary_arg},
                    stress_corr_file     = ${stress_corr_arg},
                    stress_rank_file     = ${stress_rank_arg},
                    stress_overlap_file  = ${stress_overlap_arg},
                    stress_variants_file = ${stress_variants_arg},
                    stress_loadings_file = ${stress_loadings_arg},
                    vep_transvar_file    = ${tv_arg},
                    vep_primateai_file   = ${pai_arg},
                    vep_aa2prot_file     = ${a2p_arg},
                    genomic_info_file    = ${gi_arg},
                    direction            = '${direction}'
                ),
                output_file = 'SCORING_report_${direction}.html'
            )
        "
        """
    }
}
