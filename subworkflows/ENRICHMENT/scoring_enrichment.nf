#!/usr/bin/env nextflow

/*
 * SCORING_AMI / SCORING_COMPARE
 * ────────────────────────────────────────────
 * DOMINO active-module identification (AMI) runs on the gated directional
 * 9-slice gene lists produced by SCORING_COMPUTE, with STRING used only for ID
 * mapping + per-module functional labelling. Ranked enrichment is handled by
 * SCORING_FCS_REPORT (subworkflows/ENRICHMENT/fcs.nf).
 *
 *   SCORING_AMI_REPORT     → ami_results/**, HTML (13.AMI_analysis — DOMINO modules)
 *   SCORING_COMPARE_REPORT → compare_results/**, HTML (top vs bottom: FCS only)
 */


// ─────────────────────────────────────────────────────────────────────────────
// SCORING_AMI_REPORT
// Runs 13.AMI_analysis.Rmd (DOMINO active-module identification, STRING used
// only for ID mapping + per-module functional labels) on the 9 ranked slices.
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_AMI_REPORT {
    tag "scoring_ami|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries    3

    publishDir path: "${params.outdir}/ami",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/ami/ami_results",
               mode: 'copy', overwrite: true,
               pattern: 'ami_results/**'
    publishDir path: "${params.outdir}/ami/ami_summary",
               mode: 'copy', overwrite: true,
               pattern: 'ami_summary/**'
    publishDir path: "${params.outdir}/ami/ami_plots",
               mode: 'copy', overwrite: true,
               pattern: 'ami_plots/**'
    publishDir path: "${params.outdir}/ami/ami_networks",
               mode: 'copy', overwrite: true,
               pattern: 'ami_networks/**'

    input:
    path gene_lists
    path background
    path gene_scores
    path domino_network_sif
    path domino_modules_dir
    path domino_edge_scores
    // FADE and RER each get their own DOMINO network/gene lists, built
    // against their own gene universe (see ENRICHMENT workflow) -- NO_-
    // prefixed sentinel paths when that tool didn't run this invocation.
    // DOMINO_BUILD_NETWORK/DOMINO_RUN_MODULES always emit fixed filenames
    // (network.sif / network_edge_scores.tsv / domino_modules) regardless of
    // which tool's DOMINO_MODULES call produced them -- CAAS's own 3 inputs
    // above, staged under those literal names, would otherwise collide with
    // FADE's and RER's copies the moment more than one tool actually runs in
    // the same invocation. stageAs gives FADE's and RER's copies distinct names.
    path fade_gene_lists
    // stageAs uses a non-'.txt' extension deliberately: CAAS's own gene-list
    // files are auto-detected in the Rmd via list.files(pattern="\\.txt$") at
    // the task's cwd root (minus only CAAS's own background file by name) --
    // a stray "*.txt" background file sitting at cwd root from FADE/RER would
    // otherwise be silently swept in as a phantom extra CAAS gene list.
    path fade_background, stageAs: 'fade_background.universe'
    path fade_domino_network_sif, stageAs: 'fade_network.sif'
    path fade_domino_modules_dir, stageAs: 'fade_domino_modules'
    path fade_domino_edge_scores, stageAs: 'fade_network_edge_scores.tsv'
    path rer_gene_lists
    path rer_background, stageAs: 'rer_background.universe'
    path rer_domino_network_sif, stageAs: 'rer_network.sif'
    path rer_domino_modules_dir, stageAs: 'rer_domino_modules'
    path rer_domino_edge_scores, stageAs: 'rer_network_edge_scores.tsv'
    // Explicit booleans rather than sniffing a NO_-prefixed sentinel filename:
    // stageAs above renames every FADE/RER path input to a fixed name
    // regardless of whether the real file or its sentinel arrived, so name-
    // based detection (still used for gs_arg above, whose input has no
    // stageAs) would no longer work here.
    val fade_ran
    val rer_ran

    output:
    path "13.AMI_analysis_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "ami_results/**",                     emit: ami_results,             optional: true
    path "ami_summary/**",                     emit: ami_summary,             optional: true
    path "ami_plots/**",                       emit: ami_plots,               optional: true
    path "ami_networks/**",                    emit: ami_networks,            optional: true
    path "ami_networks/module_descriptions_all_lists.tsv", emit: module_descriptions, optional: true
    path "ami_networks/term_threshold_membership.tsv",      emit: term_membership,       optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def species   = params.string_species           ?: 9606
    def net_score = params.domino_network_score_thr ?: 700
    def bg_name   = background.getName().replace("'", "\\'")
    def gs_arg    = (gene_scores.name =~ /^NO_GENE_SCORES/) ? 'NULL' : "'${gene_scores}'"

    // FADE/RER each get their own DOMINO network + gene lists (own universe,
    // see ENRICHMENT workflow) -- NO_-prefixed sentinel names mean that tool
    // didn't run this invocation, degrading every one of its params to NULL
    // so 13.AMI_analysis.Rmd skips that tool's section gracefully (same
    // pattern as gs_arg above).
    def fade_ok    = fade_ran
    def rer_ok     = rer_ran
    def fade_bg_arg      = fade_ok ? "'${fade_background}'"          : 'NULL'
    def fade_sif_arg     = fade_ok ? "'${fade_domino_network_sif}'"  : 'NULL'
    def fade_mod_arg     = fade_ok ? "'${fade_domino_modules_dir}'"  : 'NULL'
    def fade_edge_arg    = fade_ok ? "'${fade_domino_edge_scores}'"  : 'NULL'
    def fade_lists_arg   = fade_ok ? "'fade_lists'"                  : 'NULL'
    def rer_bg_arg       = rer_ok  ? "'${rer_background}'"           : 'NULL'
    def rer_sif_arg      = rer_ok  ? "'${rer_domino_network_sif}'"   : 'NULL'
    def rer_mod_arg      = rer_ok  ? "'${rer_domino_modules_dir}'"   : 'NULL'
    def rer_edge_arg     = rer_ok  ? "'${rer_domino_edge_scores}'"   : 'NULL'
    def rer_lists_arg    = rer_ok  ? "'rer_lists'"                   : 'NULL'

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

        # FADE/RER's own gene-list .txt files are staged flat by Nextflow
        # alongside CAAS's own derived *.txt files above -- move them into
        # their own subdirectories so the Rmd's per-tool sections don't pick
        # up each other's lists (each tool's list-file names are fixed/known).
        mkdir -p fade_lists rer_lists
        mv fade_top_significant.txt fade_bottom_significant.txt fade_global_significant.txt fade_lists/ 2>/dev/null || true
        mv rer_significant.txt rer_accelerating.txt rer_decelerating.txt rer_lists/ 2>/dev/null || true
    """

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                '13.AMI_analysis.Rmd',
                params = list(
                    background_file     = '${background}',
                    background_basename = '${bg_name}',
                    output_dir          = 'ami_results',
                    project_name        = '${traitname}',
                    species             = ${species},
                    domino_network_score_thr = ${net_score},
                    gene_scores_file    = ${gs_arg},
                    string_db_dir       = '${params.string_db_dir}',
                    domino_network_sif  = '${domino_network_sif}',
                    domino_modules_dir  = '${domino_modules_dir}',
                    domino_edge_scores_file = '${domino_edge_scores}',
                    fade_gene_lists_dir      = ${fade_lists_arg},
                    fade_background_file     = ${fade_bg_arg},
                    fade_domino_network_sif  = ${fade_sif_arg},
                    fade_domino_modules_dir  = ${fade_mod_arg},
                    fade_domino_edge_scores_file = ${fade_edge_arg},
                    rer_gene_lists_dir       = ${rer_lists_arg},
                    rer_background_file      = ${rer_bg_arg},
                    rer_domino_network_sif   = ${rer_sif_arg},
                    rer_domino_modules_dir   = ${rer_mod_arg},
                    rer_domino_edge_scores_file = ${rer_edge_arg}
                ),
                output_file = '13.AMI_analysis_${traitname}.html'
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
// comparison — see 13.AMI_analysis.Rmd for module-level results.
//
// Inputs (all staged flat in work dir by Nextflow):
//   fcs_all_results : fcs_results/fcs_all_results.tsv  (single file)
//
// The script sorts staged files into cmp_fcs/ before invoking
// 15.Comparison_report.Rmd.
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
    // Optional cross-angle inputs (AMI network modules, posenrich position-level
    // results) -- NO_FILE-prefixed sentinels tolerated when their module didn't run.
    path ami_module_desc,     stageAs: 'ami_module_descriptions.tsv'
    path ami_term_membership, stageAs: 'ami_term_membership.tsv'
    path posenrich_dotplot,   stageAs: 'posenrich_overall_dotplot.tsv'
    // posenrich_leading_edge.tsv (gene:position members of significant posenrich
    // terms, already position-granular) drives both the Integrated gene scorecard's
    // posenrich angle AND the Interesting Genes/Interesting Positions tables below.
    path posenrich_le,        stageAs: 'posenrich_leading_edge.tsv'
    // SCORING's own gene-level table (percentile flags + FADE/RER/Accumulation
    // significance) -- same file 12.FCS_general_report.Rmd/14.Position_enrichment_report.Rmd
    // already consume, reused here for the Interesting Genes/Positions tables.
    path fcs_stats,           stageAs: 'fcs_stats.tsv'
    // Position-level CAAS scores + FADE-site/VEP annotations -- same optional
    // inputs 14.Position_enrichment_report.Rmd's Position Characterisation section
    // already consumes, reused here for the Interesting Positions table.
    path position_scores,     stageAs: 'position_scores.tsv'
    path vep_primateai,       stageAs: 'vep_primateai.tsv'
    path vep_cosmic,          stageAs: 'vep_cosmic.tsv'
    path fade_sites_top,      stageAs: 'fade_sites_top.csv'
    path fade_sites_bottom,   stageAs: 'fade_sites_bottom.csv'
    // SCORING's published percentile slices (scoring_compute.R) -- canonical
    // gene/position 10/5/1% membership, read directly rather than re-derived
    // by this report (see 15.Comparison_report.Rmd's gene_lists_dir/
    // position_lists_dir params).
    path gene_lists,          stageAs: 'gene_lists'
    path position_lists,      stageAs: 'position_lists'

    output:
    path "15.Comparison_report_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "compare_results/**", emit: compare_results, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def fdr_thr   = params.scoring_compare_fdr   ?: params.fcs_fdr ?: 0.15
    def pperm_thr = params.fcs_pperm_thr         ?: 0.025
    def domino_module_thr = params.domino_module_thr ?: 0.05
    def top_n     = params.scoring_compare_top_n  ?: 100
    def ami_arg            = (ami_module_desc.name =~ /^NO_/)     ? 'NULL' : "'${ami_module_desc}'"
    def ami_tm_arg         = (ami_term_membership.name =~ /^NO_/) ? 'NULL' : "'${ami_term_membership}'"
    def posenrich_dot_arg  = (posenrich_dotplot.name  =~ /^NO_/) ? 'NULL' : "'${posenrich_dotplot}'"
    def posenrich_le_arg   = (posenrich_le.name  =~ /^NO_/)      ? 'NULL' : "'${posenrich_le}'"
    def fcs_stats_arg      = (fcs_stats.name =~ /^NO_/)          ? 'NULL' : "'${fcs_stats}'"
    def position_scores_arg = (position_scores.name =~ /^NO_/)   ? 'NULL' : "'${position_scores}'"
    def vep_primateai_arg  = (vep_primateai.name =~ /^NO_/)      ? 'NULL' : "'${vep_primateai}'"
    def vep_cosmic_arg     = (vep_cosmic.name =~ /^NO_/)         ? 'NULL' : "'${vep_cosmic}'"
    def fade_sites_top_arg    = (fade_sites_top.name =~ /^NO_/)    ? 'NULL' : "'${fade_sites_top}'"
    def fade_sites_bottom_arg = (fade_sites_bottom.name =~ /^NO_/) ? 'NULL' : "'${fade_sites_bottom}'"
    def gene_lists_arg     = (gene_lists.name =~ /^NO_/)     ? 'NULL' : "'${gene_lists}'"
    def position_lists_arg = (position_lists.name =~ /^NO_/) ? 'NULL' : "'${position_lists}'"

    def render_cmd = """
        Rscript -e "
            rmarkdown::render(
                '15.Comparison_report.Rmd',
                params = list(
                    fcs_dir    = 'cmp_fcs',
                    fcs_le_dir = 'cmp_fcs_le',
                    fdr_thr    = ${fdr_thr},
                    pperm_thr  = ${pperm_thr},
                    domino_module_thr = ${domino_module_thr},
                    top_n      = ${top_n},
                    traitname  = '${traitname}',
                    ami_module_desc_file      = ${ami_arg},
                    ami_term_membership_file  = ${ami_tm_arg},
                    posenrich_dotplot_file     = ${posenrich_dot_arg},
                    posenrich_leading_edge_file = ${posenrich_le_arg},
                    fcs_stats_file           = ${fcs_stats_arg},
                    position_scores_file     = ${position_scores_arg},
                    vep_primateai_file       = ${vep_primateai_arg},
                    vep_cosmic_file          = ${vep_cosmic_arg},
                    fade_sites_top_file      = ${fade_sites_top_arg},
                    fade_sites_bottom_file   = ${fade_sites_bottom_arg},
                    gene_lists_dir     = ${gene_lists_arg},
                    position_lists_dir = ${position_lists_arg}
                ),
                output_file = '15.Comparison_report_${traitname}.html'
            )
        "
    """

    def stage_cmd = """
        cp -R ${local_dir}/* .

        mkdir -p cmp_fcs cmp_fcs_le

        # Per-module FCS all-results tables (skip empty/sentinel files via -s).
        # 15.Comparison_report.Rmd derives the module from the <module>_fcs_all.tsv name.
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
