// CT Meta-CAAS Processes
// Pattern/caap_group annotation report of the discovered CAAS (the meta_caas tables themselves are written by
// CAAS_CORE_OBSERVED / CAAS_OBSERVED), plus the later join against SCORING's permulation-null significance.

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

// CAAS_SIGNIFICANCE_REPORT (note: SIGNIFICANCE, not the CT_META_CAAS process
// above -- do not confuse the two).
//
// A distinct, LATER pipeline stage than CAAS_META_CAAS_REPORT. It must
// run AFTER SCORING: it joins CT_META_CAAS's already-published
// global_meta_caas.tsv/meta_caas.tsv (pattern/caap_group breakdown, written
// upstream of SCORING in the live DAG) against SCORING's published
// position_scores.tsv (p.emp/p.adj_bh/p.adj_sam) and gene_scores.tsv
// (gene_caas_pperm/gene_caas_pperm_adj) to annotate that breakdown with the
// permulation-null significance the dead recovery_boot arm used to (badly)
// stand in for. See 16.CAAS_significance_report.Rmd for the join logic.
process CAAS_SIGNIFICANCE_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: '{significance/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path global_meta_caas
    path position_scores
    path gene_scores

    output:
    path "*.html", emit: report
    path "significance/**", emit: significance_files, optional: true
    path "significance/meta_caas_significance.tsv", emit: meta_caas_significance, optional: true
    path "significance/pattern_significance_summary.tsv", emit: pattern_significance_summary, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_META_CAAS/local"
    def outdir = "${params.outdir}/meta_caas/significance"

    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '16.CAAS_significance_report.Rmd',
                params = list(
                    global_meta_input = '${global_meta_caas}',
                    position_scores_input = '${position_scores}',
                    gene_scores_input = '${gene_scores}',
                    output_dir = '${outdir}',
                    p_emp_thr = ${params.scoring_p_emp_thr ?: 0.05},
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
                    gene_scores_input = '${gene_scores}',
                    output_dir = '${outdir}',
                    p_emp_thr = ${params.scoring_p_emp_thr ?: 0.05},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '16.CAAS_significance_report.html'
            )
        "
        """
    }
}
