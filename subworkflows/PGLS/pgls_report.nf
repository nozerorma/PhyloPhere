#!/usr/bin/env nextflow

/*
 * PGLS_REPORT
 * ───────────
 * Render an HTML summary for site-level PGLS results.
 */

process PGLS_REPORT {
    tag "pgls_report|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/characterization/pgls",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/HTML_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/characterization/pgls",
               mode: 'copy', overwrite: true,
               pattern: 'pgls_gene_summary.tsv'

    input:
    path pgls_tsv
    path site_diag_tsv
    path intval_tsv
    path extremes_tsv
    path excess_tsv

    output:
    path "PGLS_report.html", emit: report
    path "pgls_gene_summary.tsv", emit: gene_summary

    stub:
    """
    printf '%s\n' '<html><body><h1>Stub PGLS Report</h1></body></html>' > PGLS_report.html
    printf 'Gene\ttested_sites\tsignificant_sites\tbest_q_p_pgls\tbest_p_pgls\tmean_beta_top\n' > pgls_gene_summary.tsv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/PGLS/local"
    def outdir = "${params.outdir}/characterization/pgls"
    def topHits = params.pgls_top_hits_n ?: 20
    def topGenes = params.pgls_gene_top_n ?: 20
    def traitname = params.traitname ?: 'unknown_trait'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'PGLS_report.Rmd',
                params = list(
                    pgls_tsv = '${pgls_tsv}',
                    site_diag_tsv = '${site_diag_tsv}',
                    intval_tsv = '${intval_tsv}',
                    extremes_tsv = '${extremes_tsv}',
                    excess_tsv = '${excess_tsv}',
                    traitname = '${traitname}',
                    output_dir = '${outdir}',
                    top_hits_n = ${topHits},
                    top_genes_n = ${topGenes}
                ),
                output_file = 'PGLS_report.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .

        Rscript -e "
            rmarkdown::render(
                'PGLS_report.Rmd',
                params = list(
                    pgls_tsv = '${pgls_tsv}',
                    site_diag_tsv = '${site_diag_tsv}',
                    intval_tsv = '${intval_tsv}',
                    extremes_tsv = '${extremes_tsv}',
                    excess_tsv = '${excess_tsv}',
                    traitname = '${traitname}',
                    output_dir = '${outdir}',
                    top_hits_n = ${topHits},
                    top_genes_n = ${topGenes}
                ),
                output_file = 'PGLS_report.html'
            )
        "
        """
    }
}
