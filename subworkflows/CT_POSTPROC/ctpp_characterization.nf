#!/usr/bin/env nextflow
// ctpp_characterization.nf — Characterization report of the post-processed CAAS discovery.
// PhyloPhere | subworkflows/CT_POSTPROC/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_POSTPROC_REPORT: renders 8.Characterization_report.Rmd, which characterizes the
 *  discovery (genes, patterns, clusters, outliers) and the effect of the post-processing
 *  filters. Called from workflows/ct_postproc.nf.
 *
 *  Consumes:  prepared discovery, filter_summary.tsv, the directory of the published
 *             cluster files of the filter mode, gene_ensembl_file (gene lengths), gene stats
 *  Produces:  HTML report (html_reports/) and its files under postproc/ (CT_postproc_files,
 *             outliers, clusters, summary_statistics, disambiguation_characterization,
 *             postproc_inputs)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Characterization report ──────────────────────────────────────────────────


process CT_POSTPROC_REPORT {
    tag "caas_postproc_report"
    label 'process_reporting'
    publishDir path: "${params.outdir}/postproc", mode: 'copy', overwrite: true, pattern: '{CT_postproc_files/**,outliers/**,clusters/**,summary_statistics/**,disambiguation_characterization/**,postproc_inputs/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path discovery_file
    path filter_summary
    val filter_dir
    path gene_ensembl_file
    path gene_stats_file

    output:
    path "*.html", emit: report
    path "CT_postproc_files/**", emit: assets, optional: true
    path "outliers/**", emit: outliers, optional: true
    path "clusters/**", emit: clusters, optional: true
    path "summary_statistics/**", emit: summary_stats, optional: true
    path "disambiguation_characterization/**", emit: disambiguation_characterization, optional: true
    // Copy of the raw --gene_ensembl_file input. The genomic_info_file of SCORING and POSENRICH
    // is the same file (scoring.nf resolves it from params.gene_ensembl_file), so publishing it
    // here also provides the gene genomic coordinates.
    path "postproc_inputs/**", emit: gene_ensembl_input, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_POSTPROC/local"
    def discovery = discovery_file.toString()
    def filter_sum = filter_summary.toString()
    def filter_dir_path = filter_dir
    def gene_len = gene_ensembl_file.toString()
    def gene_stats = gene_stats_file.toString()
    def mode = params.caas_postproc_mode
    def outdir = "${params.outdir}/postproc"
    def extreme_thresh = params.extreme_threshold
    def iqr_mult = params.iqr_multiplier
    def gene_filter = params.gene_filter_mode
    def remove_clusters = params.remove_caas_clusters ? 'TRUE' : 'FALSE'


    if (params.use_singularity | params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        mkdir -p postproc_inputs
        [[ "\$(basename '${gene_len}')" == NO_* ]] || cp '${gene_len}' postproc_inputs/

        REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '8.Characterization_report.Rmd',
                params = list(
                    discovery_file = '${discovery}',
                    filter_summary_file = '${filter_sum}',
                    filter_dir = '${filter_dir_path}',
                    gene_ensembl_file = '${gene_len}',
                    gene_stats_file = '${gene_stats}',
                    processing_mode = '${mode}',
                    output_dir = '${outdir}',
                    extreme_threshold = ${extreme_thresh},
                    iqr_multiplier = ${iqr_mult},
                    gene_filter_mode = '${gene_filter}',
                    remove_caas_clusters = ${remove_clusters},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '8.Characterization_report.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        mkdir -p postproc_inputs
        [[ "\$(basename '${gene_len}')" == NO_* ]] || cp '${gene_len}' postproc_inputs/

        REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                '8.Characterization_report.Rmd',
                params = list(
                    discovery_file = '${discovery}',
                    filter_summary_file = '${filter_sum}',
                    filter_dir = '${filter_dir_path}',
                    gene_ensembl_file = '${gene_len}',
                    gene_stats_file = '${gene_stats}',
                    processing_mode = '${mode}',
                    output_dir = '${outdir}',
                    extreme_threshold = ${extreme_thresh},
                    iqr_multiplier = ${iqr_mult},
                    gene_filter_mode = '${gene_filter}',
                    remove_caas_clusters = ${remove_clusters},
                    seed = '${params.seed ?: 1998}'
                ),
                output_file = '8.Characterization_report.html'
            )
        "
        """
    }
}
