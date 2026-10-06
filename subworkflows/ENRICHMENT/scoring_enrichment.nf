#!/usr/bin/env nextflow
// scoring_enrichment.nf — DOMINO module report and top-versus-bottom comparison report of SCORING.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_AMI_REPORT, SCORING_COMPARE_REPORT: reports that sit on top of the SCORING
 *  gene lists and the FCS results.
 *
 *  SCORING_AMI_REPORT renders 13.AMI_analysis.Rmd: DOMINO active-module identification
 *  on the SCORING slice gene lists, with STRING used only for ID mapping and for the
 *  functional label of each module. FADE and RER enter with their own gene lists,
 *  background and DOMINO network. Ranked enrichment is in SCORING_FCS_REPORT (fcs.nf).
 *
 *  SCORING_COMPARE_REPORT renders 15.Comparison_report.Rmd: top versus bottom
 *  comparison of the FCS outputs (CAAS and RER) with the posenrich, AMI and
 *  gene/position-level evidence tables.
 *
 *  Consumes:  SCORING gene_lists/ and position_lists/, FCS and posenrich results,
 *             DOMINO networks and modules
 *  Produces:  ami/ (HTML, ami_results/, ami_summary/, ami_plots/, ami_networks/),
 *             compare/ (HTML, compare_results/)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── DOMINO module report ───────────────────────────────────────────────────────

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
    // FADE and RER each have their own DOMINO network and gene lists, built on their
    // own gene universe; NO_* sentinel paths arrive when the tool did not run.
    // DOMINO_BUILD_NETWORK and DOMINO_RUN_MODULES always write network.sif,
    // network_edge_scores.tsv and domino_modules, whichever tool they ran for. The
    // CAAS inputs above keep those names, so stageAs gives the FADE and RER copies
    // distinct names to avoid a name collision when several tools run together.
    path fade_gene_lists
    // The staged background has a '.universe' extension, not '.txt': 13.AMI_analysis.Rmd
    // detects the CAAS gene lists as every *.txt in the task root (excluding the CAAS
    // background by name), so a FADE or RER background staged as *.txt would be read
    // as an extra CAAS gene list.
    path fade_background, stageAs: 'fade_background.universe'
    path fade_domino_network_sif, stageAs: 'fade_network.sif'
    path fade_domino_modules_dir, stageAs: 'fade_domino_modules'
    path fade_domino_edge_scores, stageAs: 'fade_network_edge_scores.tsv'
    path rer_gene_lists
    path rer_background, stageAs: 'rer_background.universe'
    path rer_domino_network_sif, stageAs: 'rer_network.sif'
    path rer_domino_modules_dir, stageAs: 'rer_domino_modules'
    path rer_domino_edge_scores, stageAs: 'rer_network_edge_scores.tsv'
    // Explicit flags instead of sentinel-name detection: stageAs renames every FADE and
    // RER input, so a NO_* name never reaches the script (gene_scores has no stageAs
    // and is still detected by name, see gs_arg).
    val fade_ran
    val rer_ran

    output:
    path "13.AMI_analysis_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "ami_results/**",                     emit: ami_results,             optional: true
    path "ami_summary/**",                     emit: ami_summary,             optional: true
    path "ami_plots/**",                       emit: ami_plots,               optional: true
    path "ami_networks/**",                    emit: ami_networks,            optional: true
    path "ami_networks/ami_module_descriptions_all_tools.tsv",     emit: module_descriptions, optional: true
    path "ami_networks/ami_term_threshold_membership_all_tools.tsv", emit: term_membership,     optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def species   = params.string_species           ?: 9606
    def net_score = params.domino_network_score_thr ?: 700
    def bg_name   = background.getName().replace("'", "\\'")
    def gs_arg    = (gene_scores.name =~ /^NO_GENE_SCORES/) ? 'NULL' : "'${gene_scores}'"

    // A tool that did not run gets NULL for all its parameters, so that
    // 13.AMI_analysis.Rmd skips its section.
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

        # One gene list per slice: the Gene column of slice_<name>.tsv becomes <name>.txt
        for f in ${gene_lists}/slice_*.tsv; do
            if [ -f "\$f" ]; then
                basename=\$(basename "\$f" .tsv)
                name=\${basename#slice_}
                # First column without the header line, blank lines dropped
                tail -n +2 "\$f" | cut -f1 | { grep -v "^[[:space:]]*\$" || true; } > "\${name}.txt"
            fi
        done

        # The FADE and RER gene lists are staged flat next to the CAAS lists above;
        # moving them to their own directories keeps each tool's section of the Rmd
        # from reading the other lists (their file names are fixed).
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
                    project_name        = '${traitname}',
                    species             = ${species},
                    domino_network_score_thr = ${net_score},
                    gene_scores_file    = ${gs_arg},
                    scoring_p_emp_thr   = ${params.scoring_p_emp_thr ?: 0.05},
                    string_db_dir       = '${params.string_db_dir}',
                    string_cache_dir    = '${params.string_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/string"}',
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
                    rer_domino_edge_scores_file = ${rer_edge_arg},
                    seed                     = '${params.seed ?: 1998}'
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


// ── Top-versus-bottom comparison report ────────────────────────────────────────

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
    // One fcs_all_results.tsv per module, staged under distinct names. A module that
    // did not run arrives as the NO_FCS_ALL sentinel, which the non-empty (-s) test
    // below drops. FADE and Accumulation have no FCS ranking of their own (see fcs.nf);
    // they enter as corroboration flags on the leading edge.
    path caas_fcs,  stageAs: 'caas_fcs_all.tsv'    // CAAS (composite) FCS results
    path rer_fcs,   stageAs: 'rer_fcs_all.tsv'     // RER FCS results
    // Leading-edge tables feed the composite score (percentile concentration and
    // cross-module corroboration) that orders the tables and plots of the Cross-module
    // convergence section. NO_LEADING_EDGE sentinel when the module's report did not run.
    path caas_le,   stageAs: 'caas_leading_edge.tsv'
    path rer_le,    stageAs: 'rer_leading_edge.tsv'
    // Leading-edge composition tables exported by each module's 12.FCS_general_report.Rmd.
    path caas_le_comp, stageAs: 'caas_leading_edge_composition.tsv'
    path rer_le_comp,  stageAs: 'rer_leading_edge_composition.tsv'
    // Optional inputs from the AMI and posenrich modules; NO_* sentinels when the
    // module did not run.
    path ami_module_desc,     stageAs: 'ami_module_descriptions.tsv'
    path ami_term_membership, stageAs: 'ami_term_membership.tsv'
    path posenrich_dotplot,   stageAs: 'posenrich_overall_dotplot.tsv'
    // posenrich_leading_edge.tsv (gene:position members of the significant terms) feeds
    // the posenrich evidence of the Integrated gene scorecard and the Interesting
    // Genes and Interesting Positions tables.
    path posenrich_le,        stageAs: 'posenrich_leading_edge.tsv'
    // Leading edge table of posenrich, one row per significant term (not per cutoff),
    // exported by 14.Position_enrichment_report.Rmd for the Leading edge sub-tab of the
    // Posenrich section.
    path posenrich_le_summary, stageAs: 'posenrich_leading_edge_summary.tsv'
    // SCORING's gene-level table (percentile flags and FADE, RER and Accumulation
    // significance) for the Interesting Genes and Positions tables.
    path fcs_stats,           stageAs: 'fcs_stats.tsv'
    // Position-level CAAS scores with FADE-site and VEP annotations, for the
    // Interesting Positions table (the inputs of the Position Characterisation
    // section of 14.Position_enrichment_report.Rmd).
    path position_scores,     stageAs: 'position_scores.tsv'
    path vep_primateai,       stageAs: 'vep_primateai.tsv'
    path vep_cosmic,          stageAs: 'vep_cosmic.tsv'
    path fade_sites_top,      stageAs: 'fade_sites_top.csv'
    path fade_sites_bottom,   stageAs: 'fade_sites_bottom.csv'
    // Per-position Pfam domain/clan, UCR core/flank region and variability, and FUBAR
    // selection call (position_characterization.tsv from build_position_gmt.py): the
    // layers posenrich tests as GMTs, flattened for a Gene/Position join onto the
    // Interesting Positions table.
    path position_char,       stageAs: 'position_characterization.tsv'
    // SCORING's percentile slices (scoring_compute.R): the gene and position
    // membership of each cutoff, read by the report (gene_lists_dir and
    // position_lists_dir) instead of recomputed.
    path gene_lists,          stageAs: 'gene_lists'
    path position_lists,      stageAs: 'position_lists'
    // Inputs of the unpaired randomization null for the CAAS and RER gene-score
    // concordance: caas_perms is SCORING's caas_perms.rds (caas_corStat_byrank),
    // rer_perms is the RERconverge *.perms.rds (corStat, corRho) and gene_scores_cmp is
    // SCORING's gene_scores.tsv (observed gene_caas_score and rer_rho). NO_* sentinels
    // when unavailable.
    path caas_perms,          stageAs: 'caas_perms_cmp.rds'
    path rer_perms,           stageAs: 'rer_perms_cmp.rds'
    path gene_scores_cmp,     stageAs: 'gene_scores_cmp.tsv'

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
    def posenrich_le_summary_arg = (posenrich_le_summary.name =~ /^NO_/) ? 'NULL' : "'${posenrich_le_summary}'"
    def fcs_stats_arg      = (fcs_stats.name =~ /^NO_/)          ? 'NULL' : "'${fcs_stats}'"
    def position_scores_arg = (position_scores.name =~ /^NO_/)   ? 'NULL' : "'${position_scores}'"
    def vep_primateai_arg  = (vep_primateai.name =~ /^NO_/)      ? 'NULL' : "'${vep_primateai}'"
    def vep_cosmic_arg     = (vep_cosmic.name =~ /^NO_/)         ? 'NULL' : "'${vep_cosmic}'"
    def fade_sites_top_arg    = (fade_sites_top.name =~ /^NO_/)    ? 'NULL' : "'${fade_sites_top}'"
    def fade_sites_bottom_arg = (fade_sites_bottom.name =~ /^NO_/) ? 'NULL' : "'${fade_sites_bottom}'"
    def position_char_arg     = (position_char.name =~ /^NO_/)     ? 'NULL' : "'${position_char}'"
    def gene_lists_arg     = (gene_lists.name =~ /^NO_/)     ? 'NULL' : "'${gene_lists}'"
    def position_lists_arg = (position_lists.name =~ /^NO_/) ? 'NULL' : "'${position_lists}'"
    def caas_perms_arg     = (caas_perms.name =~ /^NO_/)      ? 'NULL' : "'${caas_perms}'"
    def rer_perms_arg      = (rer_perms.name =~ /^NO_/)       ? 'NULL' : "'${rer_perms}'"
    def gene_scores_cmp_arg = (gene_scores_cmp.name =~ /^NO_/) ? 'NULL' : "'${gene_scores_cmp}'"
    def cmp_perm_null      = (params.comparison_perm_null == null) ? true : params.comparison_perm_null
    def cmp_perm_stat      = params.comparison_perm_stat ?: 'spearman'
    def cmp_perm_topk      = params.comparison_perm_topk ?: 0.05

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
                    posenrich_leading_edge_summary_file = ${posenrich_le_summary_arg},
                    fcs_stats_file           = ${fcs_stats_arg},
                    position_scores_file     = ${position_scores_arg},
                    vep_primateai_file       = ${vep_primateai_arg},
                    vep_cosmic_file          = ${vep_cosmic_arg},
                    fade_sites_top_file      = ${fade_sites_top_arg},
                    fade_sites_bottom_file   = ${fade_sites_bottom_arg},
                    position_characterization_file = ${position_char_arg},
                    gene_lists_dir     = ${gene_lists_arg},
                    position_lists_dir = ${position_lists_arg},
                    caas_perms_file       = ${caas_perms_arg},
                    rer_perms_file        = ${rer_perms_arg},
                    gene_scores_file      = ${gene_scores_cmp_arg},
                    comparison_perm_null  = '${cmp_perm_null}',
                    comparison_perm_stat  = '${cmp_perm_stat}',
                    comparison_perm_topk  = ${cmp_perm_topk},
                    scoring_p_emp_thr  = ${params.scoring_p_emp_thr ?: 0.05},
                    seed               = '${params.seed ?: 1998}'
                ),
                output_file = '15.Comparison_report_${traitname}.html'
            )
        "
    """

    def stage_cmd = """
        cp -R ${local_dir}/* .

        mkdir -p cmp_fcs cmp_fcs_le

        # Per-module FCS tables; empty or sentinel files are skipped (-s). The Rmd
        # reads the module name from the <module>_fcs_all.tsv file name.
        for m in caas rer; do
            [ -s "\${m}_fcs_all.tsv" ] && cp "\${m}_fcs_all.tsv" "cmp_fcs/\${m}_fcs_all.tsv" || true
            [ -s "\${m}_leading_edge.tsv" ] && cp "\${m}_leading_edge.tsv" "cmp_fcs_le/\${m}_leading_edge.tsv" || true
            [ -s "\${m}_leading_edge_composition.tsv" ] && cp "\${m}_leading_edge_composition.tsv" "cmp_fcs_le/\${m}_leading_edge_composition.tsv" || true
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
