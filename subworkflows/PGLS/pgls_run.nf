#!/usr/bin/env nextflow

/*
 * SITE_PGLS
 * ─────────
 * Wrap the vendored site-level PGLS toolset and publish its TSV outputs.
 */

process SITE_PGLS {
    tag "site_pgls"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/selection/pgls",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'

    input:
    path caas_file
    path trait_file
    path tree_file

    output:
    path "site_pgls.tsv", emit: pgls_tsv
    path "site_pgls_site_diagnostics.tsv", emit: site_diag_tsv
    path "site_pgls_internal_validation.tsv", emit: intval_tsv
    path "site_pgls_extremes.tsv", emit: extremes_tsv

    stub:
    """
    printf 'Gene\tPosition\ttrait\tn_used\tn_top_aa\tn_bottom_aa\tbeta_top\tlambda_full\tlambda_null\tlrt_stat\tp_pgls\tq_p_pgls\tsig_p_pgls\n' > site_pgls.tsv
    printf 'Gene\tPosition\ttrait\tmapped_tips\tn_top_aa\tn_bottom_aa\ttrait_model\ttrait_theta\tdiag_type\tmodel\tk\ttheta\tll\taicc\tlrt_df\tlrt_stat\tlrt_p\tis_best_aicc\n' > site_pgls_site_diagnostics.tsv
    printf 'Gene\tPosition\tSubstitution\tspecies\tobserved_AA\tobserved_dir\ttrait\ttrait_value\textreme_side\tq1\tq3\n' > site_pgls_internal_validation.tsv
    printf 'species\ttrait\ttrait_value\tis_valid\tis_low\tis_high\tq1\tmedian\tq3\n' > site_pgls_extremes.tsv
    """

    script:
    def local_dir = "${baseDir}/subworkflows/PGLS/local"
    def speciesCol = params.sp_colname ?: 'species'
    def alignmentDir = params.alignment
    def alignmentFormat = params.pgls_alignment_format ?: 'auto'
    assert alignmentDir : "PGLS requires --alignment"

    def taxMapPath = params.tax_id ?: ''
    def taxMapArg = ''
    if (taxMapPath) {
        def taxMapFile = file(taxMapPath)
        assert taxMapFile.exists() : "PGLS: tax map file not found: ${taxMapPath}"
        taxMapArg = "--tax-map-file \"${taxMapFile}\""
    }

    def commonMethod = (params.discrete_method ?: 'decile').toString()
    def pglsMethod = commonMethod == 'median_sd' ? '1sd'
                   : (commonMethod in ['1sd', '2sd'] ? commonMethod : 'quantile')
    def qLower = commonMethod == 'quartile' ? '0.25'
              : commonMethod == 'quintile' ? '0.20'
              : commonMethod == 'decile' ? '0.10'
              : (params.bottom_quantile ?: '0.10')
    def qUpper = commonMethod == 'quartile' ? '0.75'
              : commonMethod == 'quintile' ? '0.80'
              : commonMethod == 'decile' ? '0.90'
              : (params.top_quantile ?: '0.90')
    def lambdaArg = params.pgls_lambda_value?.toString()?.trim() ? "--lambda-value ${params.pgls_lambda_value}" : ''
    def keepAmbiguousArg = params.pgls_keep_ambiguous ? '--keep-ambiguous' : ''
    def includeOtherArg = params.pgls_include_other ? '--include-other' : ''
    def medianSepArg = params.use_median_sep ? '--use-median-sep' : '--no-use-median-sep'

    """
    cp -R ${local_dir}/* .
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    python3 src/caas_main.py \\
        --trait-file "${trait_file}" \\
        --species-col "${speciesCol}" \\
        --trait-name "${params.traitname}" \\
        --caas-file "${caas_file}" \\
        --alignments-dir "${alignmentDir}" \\
        --alignments-format "${alignmentFormat}" \\
        --q-lower ${qLower} \\
        --q-upper ${qUpper} \\
        --extremes-method "${pglsMethod}" \\
        ${medianSepArg} \\
        ${keepAmbiguousArg} \\
        ${includeOtherArg} \\
        --tree-file "${tree_file}" \\
        --tree-format "${params.pgls_tree_format ?: 'auto'}" \\
        ${taxMapArg} \\
        --pgls-out "site_pgls.tsv" \\
        --extremes-DEBUG-out "site_pgls_extremes.tsv" \\
        --intval-out "site_pgls_internal_validation.tsv" \\
        --site-diag-out "site_pgls_site_diagnostics.tsv" \\
        --threads ${task.cpus} \\
        --site-fdr-alpha ${params.pgls_site_fdr_alpha ?: 0.05} \\
        --min-per-class ${params.pgls_min_per_class ?: 6} \\
        --phylo-select "${params.pgls_phylo_select ?: 'fixed'}" \\
        --phylo-candidates "${params.pgls_phylo_candidates ?: 'lambda'}" \\
        --phylo-bounds "${params.pgls_phylo_bounds ?: 'wide'}" \\
        --phylo-model "${params.pgls_phylo_model ?: 'lambda'}" \\
        ${lambdaArg} \\
        --ou-alpha ${params.pgls_ou_alpha ?: 0.0} \\
        --log-level "${params.pgls_log_level ?: 'INFO'}"
    """
}
