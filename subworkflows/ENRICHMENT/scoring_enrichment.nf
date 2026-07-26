#!/usr/bin/env nextflow

/*
 * SCORING_STRING / SCORING_COMPARE
 * ────────────────────────────────────────────
 * DOMINO active-module identification (AMI) runs on the gated directional
 * 9-slice gene lists produced by SCORING_COMPUTE, with STRING used only for ID
 * mapping + per-module functional labelling. Ranked enrichment is handled by
 * SCORING_FCS_REPORT (subworkflows/ENRICHMENT/fcs.nf).
 *
 *   SCORING_STRING_REPORT  → string_results/**, HTML (15.AMI_analysis — DOMINO modules)
 *   SCORING_COMPARE_REPORT → compare_results/**, HTML (top vs bottom: FCS only)
 */


// ─────────────────────────────────────────────────────────────────────────────
// SCORING_STRING_REPORT
// Runs 15.AMI_analysis.Rmd (DOMINO active-module identification, STRING used
// only for ID mapping + per-module functional labels) on the 9 ranked slices.
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_STRING_REPORT {
    tag "scoring_string|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries    3

    publishDir path: "${params.outdir}/string",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/string/string_results",
               mode: 'copy', overwrite: true,
               pattern: 'string_results/**'
    publishDir path: "${params.outdir}/string/string_summary",
               mode: 'copy', overwrite: true,
               pattern: 'string_summary/**'
    publishDir path: "${params.outdir}/string/string_plots",
               mode: 'copy', overwrite: true,
               pattern: 'string_plots/**'
    publishDir path: "${params.outdir}/string/string_networks",
               mode: 'copy', overwrite: true,
               pattern: 'string_networks/**'

    input:
    path gene_lists
    path background
    path gene_scores
    path domino_network_sif
    path domino_modules_dir

    output:
    path "15.AMI_analysis_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "string_results/**",                  emit: string_results,          optional: true
    path "string_summary/**",                  emit: string_summary,          optional: true
    path "string_plots/**",                    emit: string_plots,            optional: true
    path "string_networks/**",                 emit: string_networks,         optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def species   = params.string_species           ?: 9606
    def req_score = params.string_required_score    ?: 400
    def net_score = params.domino_network_score_thr ?: 700
    def bg_name   = background.getName().replace("'", "\\'")
    def gs_arg    = (gene_scores.name =~ /^NO_GENE_SCORES/) ? 'NULL' : "'${gene_scores}'"

    def stage_cmd = """
        cp -R ${local_dir}/* .

        # Convert slice_*.tsv to *.txt in the current directory
        for f in ${gene_lists}/slice_*.tsv; do
            if [ -f "\$f" ]; then
                basename=\$(basename "\$f" .tsv)
                name=\${basename#slice_}
                # Extract first column (Gene) except header (1st line)
                tail -n +2 "\$f" | cut -f1 | { grep -v "^[[:space:]]*\$" || true; } > "\${name}.txt"
            fi
        done
    """

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                '15.AMI_analysis.Rmd',
                params = list(
                    background_file     = '${background}',
                    background_basename = '${bg_name}',
                    output_dir          = 'string_results',
                    project_name        = '${traitname}',
                    species             = ${species},
                    required_score      = ${req_score},
                    domino_network_score_thr = ${net_score},
                    gene_scores_file    = ${gs_arg},
                    string_db_dir       = '${params.string_db_dir}',
                    domino_network_sif  = '${domino_network_sif}',
                    domino_modules_dir  = '${domino_modules_dir}'
                ),
                output_file = '15.AMI_analysis_${traitname}.html'
            )
        "
    """

    if (params.use_singularity || params.use_apptainer) {
        """
        ${stage_cmd}
        /usr/local/bin/_entrypoint.sh ${render_cmd}
        """
    } else {
        """
        ${stage_cmd}
        ${render_cmd}
        """
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SCORING_COMPARE_REPORT
// Comparative top-vs-bottom analysis across FCS outputs (CAAS + RER). DOMINO/AMI
// module output is not gene-set-ranked in the same way and is not part of this
// comparison — see 15.AMI_analysis.Rmd for module-level results.
//
// Inputs (all staged flat in work dir by Nextflow):
//   fcs_all_results : fcs_results/fcs_all_results.tsv  (single file)
//
// The script sorts staged files into cmp_fcs/ before invoking
// 12.Comparison_report.Rmd.
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_COMPARE_REPORT {
    tag "scoring_compare|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/compare",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/compare/compare_results",
               mode: 'copy', overwrite: true,
               pattern: 'compare_results/**'

    input:
    // One fcs_all_results.tsv per module — stageAs distinct names to avoid the
    // same-filename collision; absent modules arrive as the NO_FCS_ALL sentinel and
    // are filtered by the non-empty (-s) guard below + column validation in the Rmd.
    // FADE/Accumulation no longer run their own FCS ranking (unreliable — see FADE's
    // max-of-many-sites statistic and Accumulation's missing permulation null); they
    // remain available as cross-module corroboration flags on CAAS's own leading edge.
    path caas_fcs,  stageAs: 'caas_fcs_all.tsv'    // CAAS (composite) FCS results
    path rer_fcs,   stageAs: 'rer_fcs_all.tsv'     // RER FCS results
    // Leading-edge tables feed the orthogonal composite score (percentile
    // concentration + cross-module corroboration) that orders every table/plot
    // in the Cross-module convergence section. NO_LEADING_EDGE sentinel when a
    // module's report didn't run.
    path caas_le,   stageAs: 'caas_leading_edge.tsv'
    path rer_le,    stageAs: 'rer_leading_edge.tsv'

    output:
    path "12.Comparison_report_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "compare_results/**", emit: compare_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr_thr   = params.scoring_compare_fdr   ?: params.fcs_fdr ?: 0.15
    def pperm_thr = params.fcs_pperm_thr         ?: 0.025
    def top_n     = params.scoring_compare_top_n  ?: 20

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                '12.Comparison_report.Rmd',
                params = list(
                    fcs_dir    = 'cmp_fcs',
                    fcs_le_dir = 'cmp_fcs_le',
                    fdr_thr    = ${fdr_thr},
                    pperm_thr  = ${pperm_thr},
                    top_n      = ${top_n},
                    traitname  = '${traitname}'
                ),
                output_file = '12.Comparison_report_${traitname}.html'
            )
        "
    """

    def stage_cmd = """
        cp -R ${local_dir}/* .

        mkdir -p cmp_fcs cmp_fcs_le

        # Per-module FCS all-results tables (skip empty/sentinel files via -s).
        # 12.Comparison_report.Rmd derives the module from the <module>_fcs_all.tsv name.
        for m in caas rer; do
            [ -s "\${m}_fcs_all.tsv" ] && cp "\${m}_fcs_all.tsv" "cmp_fcs/\${m}_fcs_all.tsv" || true
            [ -s "\${m}_leading_edge.tsv" ] && cp "\${m}_leading_edge.tsv" "cmp_fcs_le/\${m}_leading_edge.tsv" || true
        done
    """

    if (params.use_singularity || params.use_apptainer) {
        """
        ${stage_cmd}
        /usr/local/bin/_entrypoint.sh ${render_cmd}
        """
    } else {
        """
        ${stage_cmd}
        ${render_cmd}
        """
    }
}
