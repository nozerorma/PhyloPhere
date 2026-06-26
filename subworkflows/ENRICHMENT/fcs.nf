#!/usr/bin/env nextflow

/*
 * FCS (Functional Class Scoring) report processes
 * ───────────────────────────────────────────────
 * Rank-based, threshold-free gene-set enrichment via the Wilcoxon-AUC test
 * (RERconverge::fastwilcoxGMTall) over the curated GMTs in subworkflows/ORA/dat.
 * Each process renders FCS_general.Rmd against a generic
 * stats TSV (gene + score_<ranking> + flag_<name> columns) and a universe file
 * (cleaned_background, no-signal genes floored to 0).
 *
 *   SCORING_FCS_REPORT : CAAS scoring (global/top/bottom + full cross-module flags)
 *   TOOL_FCS_REPORT    : generic, dynamic publishDir — aliased for FADE / RER
 */

// ── helper: shared param resolution ──────────────────────────────────────────
def fcs_num_g()  { params.fcs_min_genes ?: 10 }
def fcs_fdr()    { params.fcs_fdr        ?: 0.15 }   // BH (adjusted-p) cutoff
def fcs_pperm()  { params.fcs_pperm_thr  ?: 0.025 }  // permulation-p cutoff (dual threshold)
def fcs_top_n()  { params.fcs_top_n      ?: 20 }

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

    output:
    path "FCS_scoring_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "fcs_results/**",                       emit: fcs_results,     optional: true
    path "fcs_results/fcs_all_results.tsv",      emit: fcs_all_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def gmt_dir   = "${baseDir}/subworkflows/ENRICHMENT/dat"
    def num_g     = fcs_num_g()
    def fdr_thr   = fcs_fdr()
    def pperm_thr = fcs_pperm()
    def top_n     = fcs_top_n()
    def render = """
        rmarkdown::render(
            'FCS_general.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = 'Scoring_FCS_${traitname}',
                num_g         = ${num_g},
                fdr_thr       = ${fdr_thr},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${traitname}'
            ),
            output_file = 'FCS_scoring_${traitname}.html'
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
// TOOL_FCS_REPORT — generic, dynamic publishDir (aliased for FADE / RER)
// report_label drives both the publish subdir and the report filename prefix.
// ─────────────────────────────────────────────────────────────────────────────
process TOOL_FCS_REPORT {
    tag "tool_fcs|${report_label}"
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
    path annot_file

    output:
    path "${report_label}.html",            emit: report
    path "fcs_results/**",                   emit: fcs_results,     optional: true
    path "fcs_results/fcs_all_results.tsv",  emit: fcs_all_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def gmt_dir   = "${baseDir}/subworkflows/ENRICHMENT/dat"
    def num_g     = fcs_num_g()
    def fdr_thr   = fcs_fdr()
    def pperm_thr = fcs_pperm()
    def top_n     = fcs_top_n()
    def enrich_flag = (params.rer_permulation_enrichment ?: false)
    def render = """
        rmarkdown::render(
            'FCS_general.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = '${report_label}',
                num_g         = ${num_g},
                fdr_thr       = ${fdr_thr},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${params.traitname ?: "trait"}',
                enrich        = ${enrich_flag ? 'TRUE' : 'FALSE'},
                annot_file    = '${annot_file}'
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

// ─────────────────────────────────────────────────────────────────────────────
// RER_FCS_REPORT — RERconverge specific report process that takes perms_file
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
    path "${report_label}.html",            emit: report
    path "fcs_results/**",                   emit: fcs_results,     optional: true
    path "fcs_results/fcs_all_results.tsv",  emit: fcs_all_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def gmt_dir   = "${baseDir}/subworkflows/ENRICHMENT/dat"
    def num_g     = fcs_num_g()
    def fdr_thr   = fcs_fdr()
    def pperm_thr = fcs_pperm()
    def top_n     = fcs_top_n()
    def enrich_flag = (params.rer_permulation_enrichment ?: false)
    def render = """
        rmarkdown::render(
            'FCS_general.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = '${report_label}',
                num_g         = ${num_g},
                fdr_thr       = ${fdr_thr},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${params.traitname ?: "trait"}',
                perms_file    = '${perms_file}',
                enrich        = ${enrich_flag ? 'TRUE' : 'FALSE'},
                annot_file    = '${annot_file}'
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
