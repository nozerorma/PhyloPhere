#!/usr/bin/env nextflow

/*
 * FCS (Functional Class Scoring) report processes
 * ───────────────────────────────────────────────
 * Rank-based, threshold-free gene-set enrichment via the Wilcoxon-AUC test
 * (RERconverge::fastwilcoxGMTall) over the curated GMTs in subworkflows/ENRICHMENT/dat.
 * Each process renders 12.FCS_general_report.Rmd against a generic
 * stats TSV (gene + score_<ranking> + flag_<name> columns) and a universe file
 * (cleaned_background, no-signal genes floored to 0).
 *
 *   SCORING_FCS_REPORT : CAAS scoring (global/top/bottom + full cross-module flags)
 *   RER_FCS_REPORT     : RERconverge-specific report process (takes perms_file)
 *
 * FADE and Accumulation do not run their own FCS ranking: FADE's statistic is a
 * max over many sites and Accumulation has no permulation null, so neither
 * supports a reliable standalone significance test. They contribute as
 * cross-module corroboration flags on CAAS's/RER's leading edge instead, and
 * FADE additionally gets its own position-level group in posenrich.
 */

// ─────────────────────────────────────────────────────────────────────────────
// SCORING_FCS_REPORT
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_FCS_REPORT {
    tag "scoring_fcs|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/fcs",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/fcs/fcs_results",
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**'

    input:
    path fcs_stats
    path universe
    path perms_file
    path gene_lists

    output:
    path "12.FCS_scoring_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "fcs_results/**",                       emit: fcs_results,      optional: true
    path "fcs_results/fcs_all_results.tsv",      emit: fcs_all_results,  optional: true
    path "fcs_results/fcs_leading_edge.tsv",     emit: fcs_leading_edge, optional: true
    path "fcs_results/fcs_leading_edge_composition.tsv", emit: fcs_leading_edge_composition, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def gmt_dir   = params.gmt_dir
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def fdr_permsum     = params.fcs_fdr_permsum     ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    // SCORING's own published gene_lists/ -- this IS the CAAS report, so
    // score_top/score_bottom here really are CAAS's gene_caas_score_top_all/
    // bottom_all (see 12.FCS_general_report.Rmd's gene_lists_dir param doc).
    def gene_lists_arg = (gene_lists.name =~ /^NO_/) ? 'NULL' : "'${gene_lists}'"
    def render = """
        rmarkdown::render(
            '12.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = 'Scoring_FCS_${traitname}',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                fdr_wilcoxon    = ${fdr_wilcoxon},
                fdr_lachenbruch = ${fdr_lachenbruch},
                fdr_permsum     = ${fdr_permsum},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${traitname}',
                perms_file    = '${perms_file}',
                gene_lists_dir = ${gene_lists_arg},
                seed          = '${params.seed ?: 1998}'
            ),
            output_file = '12.FCS_scoring_${traitname}.html'
        )
    """
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "${render}"
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "${render}"
        """
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RER_FCS_REPORT - RERconverge specific report process that takes perms_file
// ─────────────────────────────────────────────────────────────────────────────
process RER_FCS_REPORT {
    tag "rer_fcs|${report_label}"
    label 'process_reporting'

    publishDir path: { "${params.outdir}/${subpath.toLowerCase()}" },
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: { "${params.outdir}/${subpath.toLowerCase()}/fcs_results" },
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**'

    input:
    val  subpath
    path fcs_stats
    path universe
    val  report_label
    path perms_file
    path annot_file

    output:
    path "${report_label}.html",             emit: report
    path "fcs_results/**",                   emit: fcs_results,      optional: true
    path "fcs_results/fcs_all_results.tsv",  emit: fcs_all_results,  optional: true
    path "fcs_results/fcs_leading_edge.tsv", emit: fcs_leading_edge, optional: true
    path "fcs_results/fcs_leading_edge_composition.tsv", emit: fcs_leading_edge_composition, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def gmt_dir   = params.gmt_dir
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def fdr_permsum     = params.fcs_fdr_permsum     ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    def render = """
        rmarkdown::render(
            '12.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = '${report_label}',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                fdr_wilcoxon    = ${fdr_wilcoxon},
                fdr_lachenbruch = ${fdr_lachenbruch},
                fdr_permsum     = ${fdr_permsum},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${params.traitname ?: "trait"}',
                perms_file    = '${perms_file}',
                annot_file    = '${annot_file}',
                seed          = '${params.seed ?: 1998}'
            ),
            output_file = '${report_label}.html'
        )
    """
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "${render}"
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "${render}"
        """
    }
}
