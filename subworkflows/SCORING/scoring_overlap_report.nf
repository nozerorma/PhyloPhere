#!/usr/bin/env nextflow

/*
 * SCORING_OVERLAP_REPORT
 * ──────────────────────
 * Render an HTML overlap report comparing the top-phenotype and bottom-
 * phenotype CAAS scoring results.  Calls SCORING_overlap_report.Rmd.
 *
 * Inputs
 * ──────
 *   pos_scores_top    : path — position_scores.tsv (top direction)
 *   pos_scores_bottom : path — position_scores.tsv (bottom direction)
 *   gene_scores_top   : path — gene_scores.tsv     (top direction)
 *   gene_scores_bottom: path — gene_scores.tsv     (bottom direction)
 *   genomic_info      : path — gene genomic coords TSV (NO_GENOMIC_INFO sentinel when absent)
 *
 * Outputs
 * ───────
 *   report : SCORING_report_overlap.html
 */

process SCORING_OVERLAP_REPORT {
    tag "scoring_overlap|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: { "${params.outdir}/scoring/overlap" },
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path pos_scores_top,    stageAs: 'top_position_scores.tsv'
    path pos_scores_bottom, stageAs: 'bottom_position_scores.tsv'
    path gene_scores_top,   stageAs: 'top_gene_scores.tsv'
    path gene_scores_bottom, stageAs: 'bottom_gene_scores.tsv'
    path genomic_info        // optional — NO_GENOMIC_INFO sentinel when absent

    output:
    path "SCORING_report_overlap.html", emit: report

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local"
    def outdir     = "${params.outdir}/scoring"
    def traitname  = params.traitname ?: 'unknown_trait'
    def top_pct    = params.scoring_position_top_pct  ?: 0.10
    def top5_pct   = params.scoring_position_top5_pct ?: 0.05
    def top1_pct   = params.scoring_position_top1_pct ?: 0.01
    def g_top_pct  = params.scoring_gene_top_pct      ?: 0.10
    def g_top5_pct = params.scoring_gene_top5_pct     ?: 0.05
    def g_top1_pct = params.scoring_gene_top1_pct     ?: 0.01
    def win_size   = params.scoring_window_size_bp    ?: 1000000
    def gi_arg     = (genomic_info.name =~ /^NO_GENOMIC_INFO/) ? 'NULL' : "'${genomic_info}'"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'SCORING_overlap_report.Rmd',
                params = list(
                    pos_scores_top_file    = '${pos_scores_top}',
                    pos_scores_bottom_file = '${pos_scores_bottom}',
                    gene_scores_top_file   = '${gene_scores_top}',
                    gene_scores_bottom_file= '${gene_scores_bottom}',
                    traitname              = '${traitname}',
                    top_pct                = ${top_pct},
                    top5_pct               = ${top5_pct},
                    top1_pct               = ${top1_pct},
                    gene_top_pct           = ${g_top_pct},
                    gene_top5_pct          = ${g_top5_pct},
                    gene_top1_pct          = ${g_top1_pct},
                    genomic_info_file      = ${gi_arg},
                    window_size_bp         = ${win_size}
                ),
                output_file = 'SCORING_report_overlap.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'SCORING_overlap_report.Rmd',
                params = list(
                    pos_scores_top_file    = '${pos_scores_top}',
                    pos_scores_bottom_file = '${pos_scores_bottom}',
                    gene_scores_top_file   = '${gene_scores_top}',
                    gene_scores_bottom_file= '${gene_scores_bottom}',
                    traitname              = '${traitname}',
                    top_pct                = ${top_pct},
                    top5_pct               = ${top5_pct},
                    top1_pct               = ${top1_pct},
                    gene_top_pct           = ${g_top_pct},
                    gene_top5_pct          = ${g_top5_pct},
                    gene_top1_pct          = ${g_top1_pct},
                    genomic_info_file      = ${gi_arg},
                    window_size_bp         = ${win_size}
                ),
                output_file = 'SCORING_report_overlap.html'
            )
        "
        """
    }
}
