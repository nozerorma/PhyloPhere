#!/usr/bin/env nextflow

/*
 * SITE_PGLS
 * ─────────
 * Wrap the vendored site-level PGLS toolset and publish its TSV outputs.
 */

process SITE_PGLS {
    tag "site_pgls"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/characterization/pgls",
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
    path "pgls_excess_gene_trait.tsv", emit: excess_tsv

    stub:
    """
    printf 'Gene\tPosition\ttrait\tCAAP_Group\tchange_side\tis_caas_position\tn_used\tn_top_aa\tn_bottom_aa\tbeta_top\tlambda_full\tlambda_null\tlrt_stat\tp_pgls_two_sided\tp_pgls_pos\tp_pgls_neg\tq_p_pgls_two_sided\tsig_p_pgls_two_sided\tq_p_pgls_pos\tsig_p_pgls_pos\tq_p_pgls_neg\tsig_p_pgls_neg\tp_pgls_directional\tq_p_pgls_directional\tsig_p_pgls_directional\tpgls_decile_score\tpgls_weighted_decile\tpgls_char_score\n' > site_pgls.tsv
    printf 'Gene\tPosition\ttrait\tmapped_tips\tn_top_aa\tn_bottom_aa\ttrait_model\ttrait_theta\tdiag_type\tmodel\tk\ttheta\tll\taicc\tlrt_df\tlrt_stat\tlrt_p\tis_best_aicc\n' > site_pgls_site_diagnostics.tsv
    printf 'Gene\tPosition\tSubstitution\tchange_side\tCAAP_Group\tis_caas_position\tspecies\tobserved_AA_raw\tobserved_AA\tobserved_dir\ttrait\ttrait_value\textreme_side\tq1\tq3\n' > site_pgls_internal_validation.tsv
    printf 'species\ttrait\ttrait_value\tis_valid\tis_low\tis_high\tq1\tmedian\tq3\n' > site_pgls_extremes.tsv
    printf 'trait\tGene\tCAAP_Group\tN_testable\tK_sig_testable\tn_testable_caas\tk_sig_caas\trate_testable\trate_caas\texcess_ratio\tp_hyper_excess\tq_hyper_excess\tsig_hyper_excess\n' > pgls_excess_gene_trait.tsv
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

    def caasTraitPath = params.pgls_caas_trait_file?.trim()
                        ?: "${params.outdir}/data_exploration/2.CT/1.Traitfiles/traitfile.tab"
    def caasTraitFile = file(caasTraitPath)
    def caasTraitArg  = caasTraitFile.exists() ? "--caas-trait-file \"${caasTraitFile}\"" : ''

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
    def toyMode = params.toy_mode ? 'true' : 'false'
    def toyN    = (params.toy_n ?: 50) as int
    // pgls_enrichment=true forces all-position mode so the hypergeometric
    // excess test has a non-CAAS background; otherwise honour pgls_position_mode.
    def positionMode = params.pgls_enrichment ? 'all_positions_for_caas_genes'
                     : (params.pgls_position_mode ?: 'caas_only')

    """
    cp -R ${local_dir}/* .
    mkdir -p src/biochem
    cp ${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/biochem/grouping.py src/biochem/grouping.py
    find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name '*.pyc' -delete 2>/dev/null || true

    _ali_dir="${alignmentDir}"

    # Toy mode: subsample alignments, ensuring all required genes are included.
    if [[ "${toyMode}" == "true" ]]; then
        mkdir -p _ali_toy
        
        # Step 1: Link alignments matching genes in CAAS file
        grep -v "^gene" "${caas_file}" | cut -f1 | sort -u > _required_genes.txt
        
        for _gene in \$(cat _required_genes.txt); do
            for _file in "\$_ali_dir"/\${_gene}.*; do
                [[ -f "\$_file" ]] && ln -sf "\$_file" "_ali_toy/\$(basename "\$_file")"
            done
        done
        
        # Step 2: Count current symlinks and add random ones to reach 150
        for _file in "\$_ali_dir"/*; do
            _current=\$(ls -1 _ali_toy 2>/dev/null | wc -l)
            [[ \$_current -ge 150 ]] && break
            _basefile=\$(basename "\$_file")
            [[ ! -e "_ali_toy/\$_basefile" ]] && ln -sf "\$_file" "_ali_toy/\$_basefile"
        done
        
        _ali_dir="\$(pwd)/_ali_toy"
    fi

    python3 src/caas_main.py \\
        --trait-file "${trait_file}" \\
        --species-col "${speciesCol}" \\
        --trait-name "${params.traitname}" \\
        --caas-file "${caas_file}" \\
        --alignments-dir "\$_ali_dir" \\
        --alignments-format "${alignmentFormat}" \\
        --q-lower ${qLower} \\
        --q-upper ${qUpper} \\
        --extremes-method "${pglsMethod}" \\
        ${medianSepArg} \\
        ${keepAmbiguousArg} \\
        ${includeOtherArg} \\
        --position-mode "${positionMode}" \\
        --tree-file "${tree_file}" \\
        --tree-format "${params.pgls_tree_format ?: 'auto'}" \\
        ${taxMapArg} \\
        ${caasTraitArg} \\
        ${params.n_trait ? "--n-trait \"${params.n_trait}\"" : ''} \\
        --pgls-out "site_pgls.tsv" \\
        --extremes-DEBUG-out "site_pgls_extremes.tsv" \\
        --intval-out "site_pgls_internal_validation.tsv" \\
        --site-diag-out "site_pgls_site_diagnostics.tsv" \\
        --excess-out "pgls_excess_gene_trait.tsv" \\
        --threads ${task.cpus} \\
        --site-fdr-alpha ${params.pgls_site_fdr_alpha ?: 0.05} \\
        --min-per-class ${params.pgls_min_per_class ?: 4} \\
        --phylo-select "${params.pgls_phylo_select ?: 'fixed'}" \\
        --phylo-candidates "${params.pgls_phylo_candidates ?: 'lambda'}" \\
        --phylo-bounds "${params.pgls_phylo_bounds ?: 'wide'}" \\
        --phylo-model "${params.pgls_phylo_model ?: 'lambda'}" \\
        ${lambdaArg} \\
        --ou-alpha ${params.pgls_ou_alpha ?: 0.0} \\
        --log-level "${params.pgls_log_level ?: 'INFO'}"
    """
}
