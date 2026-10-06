#!/usr/bin/env nextflow
// scoring_report.nf — Render the HTML report of the position-level and gene-level CAAS scores.
// PhyloPhere | subworkflows/SCORING/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_REPORT: renders 11.Scoring_report.Rmd for one trait from the tables of
 *  SCORING_COMPUTE and the optional evidence files.
 *
 *  Consumes:  position_scores.tsv, gene_scores.tsv, gene_correlations.tsv; optional
 *             per-site FADE tables (top and bottom), gene coordinates, the CAAS
 *             permulation null (caas_perms.rds), filtered_discovery.tsv. The
 *             optional inputs are NO_* sentinel files when absent and reach the report
 *             as NULL; caas_pos_sample, caas_pos_quantiles and background_file are
 *             passed to the report but its body does not read them.
 *  Produces:  11.Scoring_report_<trait>.html, published to scoring/ and html_reports/
 *             (the report also writes position_scores_annotated.tsv and the genomic
 *             window tables into scoring/)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Scoring report ─────────────────────────────────────────────────────────────

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
    path fade_site_top_file  // optional: per-site FADE table, top direction (NO_FADE_SITE_TOP sentinel when absent)
    path fade_site_bot_file  // optional: per-site FADE table, bottom direction (NO_FADE_SITE_BOT sentinel when absent)
    path genomic_info        // optional: gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)
    path caas_perms          // optional: caas_perms.rds, the ASR and CAAS permulation null (NO_CAAS_PERMS or NO_FILE sentinel when absent)
    path caas_pos_sample     // optional: per-scheme null sample (NO_CAAS_POS_SAMPLE sentinel); passed to the report, which does not read it
    path caas_pos_quantiles  // optional: per (cycle, scheme) null quantiles (NO_CAAS_POS_QUANTILES sentinel); passed to the report, which does not read it
    path filtered_discovery  // filtered_discovery.tsv: observed asr_path_score per (gene, position, scheme), read by the report for the biochemistry and convergence-type panels
    path background_file     // background gene list (NO_BACKGROUND or NO_FILE sentinel when absent); passed to the report, which does not read it

    output:
    path "11.Scoring_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def local_dir      = "${baseDir}/subworkflows/SCORING/local"
    def outdir         = "${params.outdir}/scoring"
    def traitname      = params.traitname ?: 'unknown_trait'
    def top_pct        = params.scoring_position_top_pct   ?: 0.10
    def g_top_pct      = params.scoring_gene_top_pct       ?: 0.10
    // Optional inputs: a sentinel file becomes the R literal NULL, a real file its quoted name
    def fs_top_arg = (fade_site_top_file.name =~ /^NO_FADE_SITE_TOP/) ? 'NULL' : "'${fade_site_top_file}'"
    def fs_bot_arg = (fade_site_bot_file.name =~ /^NO_FADE_SITE_BOT/) ? 'NULL' : "'${fade_site_bot_file}'"

    def gi_arg  = (genomic_info.name  =~ /^NO_GENOMIC_INFO/)  ? 'NULL' : "'${genomic_info}'"
    def perms_arg = (caas_perms.name =~ /^NO_CAAS_PERMS|^NO_FILE/) ? 'NULL' : "'${caas_perms}'"
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
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    genomic_info_file    = ${gi_arg},
                    caas_perms_file      = ${perms_arg},
                    caas_pos_sample_file = ${pos_sample_arg},
                    caas_pos_quantiles_file = ${pos_quantiles_arg},
                    filtered_discovery_file = ${filt_disc_arg},
                    background_file      = ${bg_file_arg},
                    scoring_p_emp_thr    = ${params.scoring_p_emp_thr ?: 0.05},
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
                    fade_site_top_file   = ${fs_top_arg},
                    fade_site_bot_file   = ${fs_bot_arg},
                    genomic_info_file    = ${gi_arg},
                    caas_perms_file      = ${perms_arg},
                    caas_pos_sample_file = ${pos_sample_arg},
                    caas_pos_quantiles_file = ${pos_quantiles_arg},
                    filtered_discovery_file = ${filt_disc_arg},
                    background_file      = ${bg_file_arg},
                    scoring_p_emp_thr    = ${params.scoring_p_emp_thr ?: 0.05},
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
