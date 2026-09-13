// CT Signification Processes
// Perform significance testing via hypergeometric and permutation analysis

process CAAS_SIGNIFICATION_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/signification", mode: 'copy', overwrite: true, pattern: '{CAAS_signification_files/**,meta_caas/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path discovery_input
    path background_input

    output:
    path "*.html", emit: report
    path "CAAS_signification_files/**", emit: assets, optional: true
    path "meta_caas/**", emit: meta_caas, optional: true
    path "meta_caas/global_meta_caas.tsv", emit: global_meta_caas, optional: true

    script:
    def caap_mode_r = params.caap_mode ? 'TRUE' : 'FALSE'
    def local_dir = "${baseDir}/subworkflows/CT_SIGNIFICATION/local"
    def outdir = "${params.outdir}/signification"


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '7.CT_signification.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '7.CT_signification.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        # Render R Markdown report
        Rscript -e "
            rmarkdown::render(
                '7.CT_signification.Rmd',
                params = list(
                    discovery_input = '${discovery_input}',
                    background_input = '${background_input}',
                    output_dir = '${outdir}',
                    caap_mode = ${caap_mode_r},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '7.CT_signification.html'
            )
        "
        """
    }
}

// CAAS_SIGNIFICANCE_REPORT (note: SIGNIFICANCE, not SIGNIFICATION -- do not
// confuse with CAAS_SIGNIFICATION_REPORT above).
//
// A distinct, LATER pipeline stage than CAAS_SIGNIFICATION_REPORT. It must
// run AFTER SCORING: it joins CT_SIGNIFICATION's already-published
// global_meta_caas.tsv/meta_caas.tsv (pattern/caap_group breakdown, written
// upstream of SCORING in the live DAG) against SCORING's published
// position_scores.tsv (p.emp/p.emp_adj) and gene_scores.tsv
// (gene_caas_pperm/gene_caas_pperm_adj) to annotate that breakdown with the
// permulation-null significance the dead recovery_boot arm used to (badly)
// stand in for. See 16.CAAS_significance_report.Rmd for the join logic.
process CAAS_SIGNIFICANCE_REPORT {
    label 'process_reporting'
    publishDir path: "${params.outdir}/signification", mode: 'copy', overwrite: true, pattern: '{significance/**}'
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
    def local_dir = "${baseDir}/subworkflows/CT_SIGNIFICATION/local"
    def outdir = "${params.outdir}/signification/significance"

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
                    p_emp_thr = ${params.scoring_p_emp_thr ?: 0.1},
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
                    p_emp_thr = ${params.scoring_p_emp_thr ?: 0.1},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '16.CAAS_significance_report.html'
            )
        "
        """
    }
}
