// ctpp_meta_caas.nf — Reports on the pattern/caap_group annotation of the CAAS and on its permulation significance.
// PhyloPhere | subworkflows/CT_META_CAAS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_META_CAAS: two R Markdown reports. The meta_caas tables themselves are written by
 *  CAAS_CORE_OBSERVED / CAAS_OBSERVED, not here.
 *
 *  Consumes:  CAAS_META_CAAS_REPORT: discovery file and background gene list
 *             CAAS_SIGNIFICANCE_REPORT: per-position CAAS table (global_meta_caas.tsv or the
 *             post-processed discovery), position_scores.tsv and gene_scores.tsv of SCORING
 *  Produces:  HTML reports (html_reports/), meta_caas/CAAS_pattern_annotation_files/ and
 *             meta_caas/significance/ (meta_caas_significance.tsv, pattern_significance_summary.tsv)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Pattern annotation report ────────────────────────────────────────────────

// Renders 7.CAAS_pattern_annotation.Rmd; called from workflows/ct_meta_caas.nf.
process CAAS_META_CAAS_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: 'CAAS_pattern_annotation_files/**'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path discovery_input
    path background_input

    output:
    path "*.html", emit: report
    path "CAAS_pattern_annotation_files/**", emit: assets, optional: true

    script:
    def caap_mode_r = params.caap_mode ? 'TRUE' : 'FALSE'
    def local_dir = "${baseDir}/subworkflows/CT_META_CAAS/local"
    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '7.CAAS_pattern_annotation.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = '7.CAAS_pattern_annotation.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        Rscript -e "
            rmarkdown::render(
                '7.CAAS_pattern_annotation.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    caap_mode = ${caap_mode_r}
                ),
                output_file = '7.CAAS_pattern_annotation.html'
            )
        "
        """
    }
}


// ── Significance report ──────────────────────────────────────────────────────

// Renders 16.CAAS_significance_report.Rmd; called from main.nf. It runs after SCORING, because it joins the
// p.emp, p.adj_bh, p.emp_fact and p.adj_bh_fact of position_scores.tsv onto the CAAS table; the join logic is in
// the Rmd. The p-values are informative and no threshold is applied. Not to be confused with
// CAAS_META_CAAS_REPORT, which runs upstream of SCORING.
process CAAS_SIGNIFICANCE_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: '{significance/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path global_meta_caas
    path position_scores

    output:
    path "*.html", emit: report
    path "significance/**", emit: significance_files, optional: true
    path "significance/meta_caas_significance.tsv", emit: meta_caas_significance, optional: true
    path "significance/pattern_significance_summary.tsv", emit: pattern_significance_summary, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_META_CAAS/local"
    def outdir = "${params.outdir}/meta_caas/significance"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '16.CAAS_significance_report.Rmd',
                params = list(
                    global_meta_input = '${global_meta_caas}',
                    position_scores_input = '${position_scores}',
                    output_dir = '${outdir}',
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '16.CAAS_significance_report.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        Rscript -e "
            rmarkdown::render(
                '16.CAAS_significance_report.Rmd',
                params = list(
                    global_meta_input = '${global_meta_caas}',
                    position_scores_input = '${position_scores}',
                    output_dir = '${outdir}',
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '16.CAAS_significance_report.html'
            )
        "
        """
    }
}
