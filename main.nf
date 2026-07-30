#!/usr/bin/env nextflow

/*
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#                                                                                    
#                                      
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: main.nf
#
*/

/*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Unlock the secrets of evolutionary relationships with Phylophere! 🌳🔍 This Nextflow pipeline
* packs a powerful punch, offering a comprehensive suite of phylogenetic comparative tools
* and analyses. Dive into the world of evolutionary biology like never before and elevate
* your research to new heights! 🚀🧬 #Phylophere #EvolutionaryInsights #NextflowPipeline
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

version = "2.0.0"

// Display input parameters
log.info """

PHYLOPHERE - NF PIPELINE  ~  version ${version}
=============================================

PHYLOPHERE: A Nextflow pipeline including a complete set
of phylogenetic comparative tools and analyses for Phenome-Genome studies

Author:         Miguel Ramon (miguel.ramon@upf.edu)


"""

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  NAMED WORKFLOW FOR PIPELINE: This section includes the main workflows.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include {HELP} from './workflows/help.nf'
include {CT} from './workflows/ct.nf'
include {RER_MAIN} from './workflows/rerconverge.nf'
include {REPORTING} from './workflows/reporting.nf'
include {CONTRAST_SELECTION} from './workflows/contrast_selection.nf'
include {CT_SIGNIFICATION} from './workflows/ct_signification.nf'
include {CT_POSTPROC} from './workflows/ct_postproc.nf'
include {CT_DISAMBIGUATION} from './workflows/ct_disambiguation.nf'
include {CT_ACCUMULATION} from './workflows/ct_accumulation.nf'
include {FADE}           from './workflows/fade.nf'
include {FADE_REPORT as FADE_REPORT_PRECOMP_TOP; FADE_REPORT as FADE_REPORT_PRECOMP_BOTTOM} from './subworkflows/FADE/fade_report.nf'
include {FADE_GENE_LISTS as FADE_GENE_LISTS_PRECOMP_TOP; FADE_GENE_LISTS as FADE_GENE_LISTS_PRECOMP_BOTTOM} from './subworkflows/FADE/fade_gene_lists.nf'
include {SELECTION_PREP} from './subworkflows/SELECTION/selection_prep.nf'
include {VEP}            from './workflows/vep.nf'
include {SCORING}        from './workflows/scoring.nf'
include {CAAS_PERMULATION} from './subworkflows/CT/caas_permulation.nf'
include {ENRICHMENT}      from './workflows/enrichment.nf'

// Workflow-map helper logic lives in lib/WorkflowMap.groovy (auto-loaded by Nextflow)

// Post-completion: write workflow_map.html again after all publishDir copies are done.
workflow.onComplete {
    try {
        // Resolve to an absolute canonical path so the file is always written
        // to the correct location regardless of JVM working directory at hook time.
        def outdirRaw = params.outdir ? params.outdir.toString() : "${workflow.projectDir}/out"
        def outdirAbs = new File(outdirRaw).canonicalPath
        def ctx = WorkflowMap.buildCtx(outdirAbs, params, workflow)
        def outdirFile = new File(outdirAbs)
        if (!outdirFile.exists()) outdirFile.mkdirs()
        def html = WorkflowMap.buildWorkflowMapHtml(ctx)

        // Primary artifact name
        def mapTarget = new File(outdirFile, 'workflow_map.html')
        log.info "[FINAL_HTML] Workflow map target: ${mapTarget.absolutePath}"
        mapTarget.text = html

        // Explicit completion marker so users can quickly verify final HTML generation.
        def markerTarget = new File(outdirFile, 'workflow_html.done')
        markerTarget.text = """status=ok
workflow_map=${mapTarget.absolutePath}
generated_at=${new Date().format("yyyy-MM-dd'T'HH:mm:ssXXX")}
"""

        log.info "[FINAL_HTML] Workflow map generated: ${mapTarget.absolutePath}"
        log.info "[FINAL_HTML] Completion marker generated: ${markerTarget.absolutePath}"
    } catch (Throwable t) {
        log.warn "Could not generate final workflow map HTML: [${t.class.simpleName}] ${t.message ?: '(null message)'}"
        t.printStackTrace()
    }
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RUN PHYLOPHERE ANALYSIS: This section initiates the main Phylophere workflow.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow {

    // Check if --help is provided
    if (params.help) {
        HELP ()
    } else {
        // Run any combination of tools requested
        def ran_any = false
        def reporting_results = null
        def contrast_out = null
        // Populated by the precomputed-FADE-input block below (json_dir -> summary/site
        // TSVs) so SCORING can reuse them without a separate --scoring_fade_* path.
        def fade_precomp_top_out = null
        def fade_precomp_bot_out = null

        if (params.reporting && !params.contrast_selection) {
            reporting_results = REPORTING()
            ran_any = true
        }
        def ct_results
        if (params.ct_tool) {
            if (params.contrast_selection) {
                contrast_out = CONTRAST_SELECTION()

                // Hard stop: if CHECK_MIN_CONTRASTS emits low_contrasts.skip,
                // terminate the current trait run gracefully (exit 0).
                contrast_out.low_contrasts_skip.view { skip_file ->
                    exit 0, "Minimum contrast threshold not met for trait '${params.traitname ?: 'unknown'}' (flag: ${skip_file}). Stopping pipeline gracefully."
                }

                ct_results = CT(contrast_out.trait_file_out, contrast_out.bootstrap_trait_file_out, contrast_out.tree_file_out)
            } else {
                def trait_file_in = null
                def bootstrap_trait_file_in = null
                def tree_file_in = null
                ct_results = CT (trait_file_in, bootstrap_trait_file_in, tree_file_in)
            }
            ran_any = true
        }
        if (params.contrast_selection && !params.ct_tool) {
            contrast_out = CONTRAST_SELECTION()
            ran_any = true
        }
        def signification_results = null
        def disambiguation_results = null
        def postproc_results = null

        // Track which sub-tools actually ran so we can pass null (not Channel.empty())
        // to downstream workflows when a tool didn't produce output.
        // Channel.empty() is truthy in Groovy, so if() guards inside sub-workflows
        // would take the "use CT output" branch but the channel would never emit,
        // causing all downstream .ifEmpty{} fallbacks to fire incorrectly.
        // Guard: params.ct_tool may be Boolean true if --ct_tool "" was passed
        // on some shells; use instanceof check before calling .split().
        def ct_tools_ran      = (params.ct_tool instanceof String && params.ct_tool)
                                    ? params.ct_tool.split(',').collect { it.trim() } : []
        def ran_discovery     = ct_tools_ran.contains('discovery')
        def ran_bootstrap     = ct_tools_ran.contains('bootstrap')

        // Signification always runs downstream of CAAStools when bootstrap produced
        // permutation output (or a standalone --bootstrap_from directory is supplied).
        // No separate --ct_signification toggle: it is implied by running bootstrap.
        def run_signification = ran_bootstrap || params.bootstrap_from

        // Stable channel references for CT_POSTPROC outputs used by multiple consumers.
        // Populated inside the ct_postproc block when --ct_postproc is enabled.
        def pp_cleaned_bg     = null   // cleaned_background_main (single file, value channel)

        if (run_signification) {
            // Only pass CT channels when the corresponding tool actually ran.
            // Pass null (not Channel.empty()) when absent so the if(channel) guard
            // inside CT_SIGNIFICATION correctly detects absence and falls back to params.
            def discovery_ch        = (ct_results && ran_discovery) ? ct_results.discovery_file   : null
            def background_genes_ch = (ct_results && ran_discovery) ? ct_results.background_genes  : null
            def bootstrap_ch        = (ct_results && ran_bootstrap) ? ct_results.bootstrap_file    : null

            signification_results = CT_SIGNIFICATION(discovery_ch, background_genes_ch, bootstrap_ch)
            ran_any = true
        }

        if (params.ct_disambiguation) {
            if (!run_signification && !params.signification_from) {
                error "CT disambiguation requires upstream signification (run CAAStools bootstrap via --ct_tool ...,bootstrap) or a standalone metadata file (--signification_from)."
            }

            // Forward both possible signification metadata artifacts; CT_DISAMBIGUATION
            // will prefer global_meta_caas.tsv when present and otherwise accept
            // the per-run meta_caas.tsv fallback.
            def meta_for_disambiguation = signification_results
                ? signification_results.signification_global_meta.mix(signification_results.signification_meta_caas)
                : null
            def trait_for_disambiguation = ct_results ? ct_results.trait_file : Channel.empty()
            def tree_for_disambiguation = ct_results ? ct_results.tree_file : Channel.empty()

            disambiguation_results = CT_DISAMBIGUATION(meta_for_disambiguation, trait_for_disambiguation, tree_for_disambiguation)
            ran_any = true
        }


        if (params.ct_postproc) {
            if (!params.ct_disambiguation && !params.disambiguation_input) {
                error "CT post-processing now runs downstream of disambiguation. Enable --ct_disambiguation or provide --disambiguation_input (caas_convergence_master.csv)."
            }

            // Post-processing is downstream from disambiguation; consume disambiguation master CSV when available
            // Pass null (not Channel.empty()) when there is no upstream result so that the
            // if(channel) guard inside CT_POSTPROC correctly detects absence and falls back
            // to --disambiguation_input / --background_input params (same pattern as CT_SIGNIFICATION).
            def disambiguation_ch = disambiguation_results ? disambiguation_results.master_csv : null
            // Only wire raw background channels when discovery actually ran; otherwise pass
            // null/Channel.empty() so CT_POSTPROC falls back to --background_input param.
            def background_ch       = (ct_results && ran_discovery) ? ct_results.background_file_raw : Channel.empty()
            def background_genes_ch = (ct_results && ran_discovery) ? ct_results.background_genes    : null
            def bootstrap_ch = Channel.empty() // retained for CT_POSTPROC signature compatibility
            // Pass full ct_disambiguation/ directory for ASR robustness diagnostics (null = standalone mode)
            def disambiguation_dir_ch = disambiguation_results ? disambiguation_results.results_dir : null

            postproc_results = CT_POSTPROC(disambiguation_ch, background_ch, background_genes_ch, bootstrap_ch, disambiguation_dir_ch)
            ran_any = true

            // Capture postproc outputs as reusable references.
            // cleaned_background is already a value channel (single file from CAAS_BACKGROUND_CLEANUP).
            pp_cleaned_bg = postproc_results.cleaned_background
        }

        def accum_results = null

        if (params.ct_accumulation) {
            if (!params.ct_postproc && !params.accumulation_background_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone background file (--accumulation_background_input)."
            }
            if (!params.ct_postproc && !params.accumulation_caas_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone CAAS file (--accumulation_caas_input)."
            }

            // Use filtered_discovery.tsv from postproc (gene_filtering stage)
            def acc_caas_ch       = postproc_results ? postproc_results.filtered_discovery : Channel.empty()
            def acc_background_ch = pp_cleaned_bg    ?: Channel.empty()
            def acc_trait_file_ch = ct_results       ? ct_results.trait_file               : Channel.empty()

            accum_results = CT_ACCUMULATION(acc_caas_ch, acc_background_ch, acc_trait_file_ch)
            ran_any = true

        }

        if (params.vep) {
            // Pass null (not Channel.empty()) when there is no upstream postproc
            // output so VEP can correctly fall back to --vep_caas_input.
            def vep_caas_ch = postproc_results ? postproc_results.filtered_discovery : null
            VEP(vep_caas_ch)
            ran_any = true
        }

        if (params.fade) {
            // Resolve upstream channel sources for SELECTION_PREP.
            // These channels are now consumed by a single SELECTION_PREP call
            // (instead of being split separately for FADE).

            def stats_source_ch = contrast_out
                ? contrast_out.stats_file_out
                : (reporting_results ? reporting_results.stats_file : Channel.empty())

            // Tree: prefer phenotype-pruned tree from reporting; fall back to
            // CT-pruned tree from contrast_selection; otherwise SELECTION_PREP
            // falls back to params.tree via its .ifEmpty {} guard.
            def tree_source_ch = reporting_results
                ? reporting_results.pruned_tree_file
                : (contrast_out ? contrast_out.tree_file_out : Channel.empty())

            // CT discovery output for toy_mode gene reuse (null when CT didn't run)
            def ct_discovery_source_ch = (ct_results && ran_discovery)
                ? ct_results.discovery_file
                : Channel.empty()

            // Postproc directional gene lists for gene_set mode.
            // SELECTION_PREP uses all_top.txt (top+both CAAS hits) to restrict
            // the top FADE run, and all_bottom.txt for the bottom run.
            // collectFile() merges multi-phenotype emissions into a single file so
            // COLLECT_GENE_SETS inside SELECTION_PREP runs exactly once.
            def sel_pp_top_ch    = Channel.empty()
            def sel_pp_bottom_ch = Channel.empty()
            if (postproc_results) {
                def pp_gene_lists_val = postproc_results.enrichment_gene_lists_files.collect()
                sel_pp_top_ch    = pp_gene_lists_val
                    .flatMap { files -> files.findAll { f -> f.name == 'all_top.txt' } }
                    .collectFile(name: 'merged_pp_top.txt')
                sel_pp_bottom_ch = pp_gene_lists_val
                    .flatMap { files -> files.findAll { f -> f.name == 'all_bottom.txt' } }
                    .collectFile(name: 'merged_pp_bottom.txt')
            }

            // Run alignment prep ONCE for FADE.
            // SELECTION_PREP outputs value channels for species/tree files
            // (can be safely consumed by multiple downstream operators) and
            // queue channels for the per-gene filtered FASTAs.
            SELECTION_PREP(
                stats_source_ch,
                tree_source_ch,
                sel_pp_top_ch,
                sel_pp_bottom_ch,
                ct_discovery_source_ch
            )

            // Fork each per-gene fasta channel into independent copies for
            // FADE  so they don't compete on the same queue channel.
            def prep_top_fanout    = SELECTION_PREP.out.filtered_fasta_top_ch
                .multiMap { it -> fade: it}
            def prep_bottom_fanout = SELECTION_PREP.out.filtered_fasta_bottom_ch
                .multiMap { it -> fade: it}

            if (params.fade) {
                FADE(
                    prep_top_fanout.fade,
                    prep_bottom_fanout.fade,
                    SELECTION_PREP.out.tree_ch,
                    SELECTION_PREP.out.top_species,
                    SELECTION_PREP.out.bottom_species
                )
                ran_any = true
            }
        }

        // Precomputed-FADE-input path: render 6.FADE_report_{top,bottom}.html
        // directly from prior *.FADE.json results when --fade itself doesn't
        // run this invocation (mirrors rer_continuous_file's role for RER,
        // just below). FADE_REPORT already emits summary_tsv/site_tsv as a
        // side effect of rendering (optional, same process SCORING would use
        // if --fade ran live) — captured here into fade_precomp_{top,bot}_out
        // so the SCORING block further down can wire them in directly instead
        // of requiring a separately-supplied --scoring_fade_summary_*/--scoring_fade_site_*
        // TSV. Those manual params remain as a fallback inside scoring.nf for
        // the case where the caller already has the TSV but not the raw JSONs.
        if (!params.fade && (params.fade_json_dir_top || params.fade_json_dir_bottom)) {
            def precomp_top_jsons = params.fade_json_dir_top
                ? Channel.fromPath("${params.fade_json_dir_top}/*.FADE.json").collect().ifEmpty([])
                : Channel.value([])
            def precomp_bottom_jsons = params.fade_json_dir_bottom
                ? Channel.fromPath("${params.fade_json_dir_bottom}/*.FADE.json").collect().ifEmpty([])
                : Channel.value([])
            fade_precomp_top_out = FADE_REPORT_PRECOMP_TOP(Channel.value('top'), precomp_top_jsons, file('NO_FG_LIST'))
            fade_precomp_bot_out = FADE_REPORT_PRECOMP_BOTTOM(Channel.value('bottom'), precomp_bottom_jsons, file('NO_FG_LIST'))
            ran_any = true
        }

        // rer_continuous_file lets RER_MAIN render 5.RERconverge_report.html from a
        // precomputed *.continuous.output RDS even when --rer_tool itself is off
        // this invocation (see rerconverge.nf's own params.rer_continuous_file branch).
        if (params.rer_tool || params.rer_continuous_file) {
            // NOTE: RER_TRAIT requires the original phenotype file (with proper column
            // headers), NOT the caastools traitfile (headerless 3-col format).
            // Always pass Channel.empty() so RER_MAIN falls back to --my_traits.
            def rer_traitfile_ch = Channel.empty()
            RER_MAIN(
                rer_traitfile_ch,
                Channel.empty(),
                Channel.empty()
            )
            ran_any = true
        }

        println "DEBUG: params.traitname = '${params.traitname}'"
        if (params.scoring) {
            if (!params.ct_postproc && !params.scoring_postproc_input) {
                error "SCORING requires CT post-processing output (--ct_postproc) or --scoring_postproc_input."
            }

            // Wire upstream outputs into SCORING. Pass null (not Channel.empty())
            // when a module didn't run so if(channel) guards detect absence correctly.
            def scoring_postproc_ch      = postproc_results ? postproc_results.filtered_discovery : null
            // Precomputed FADE JSONs (--fade_json_dir_top/_bottom, no --fade this run)
            // feed the exact same summary_tsv/site_tsv SCORING needs, via
            // fade_precomp_{top,bot}_out captured above — so a --scoring run against a
            // prior --fade run's JSON output no longer needs separately pre-built
            // --scoring_fade_summary_*/--scoring_fade_site_* TSVs.
            def scoring_fade_top_ch      = params.fade ? FADE.out.summary_top    : fade_precomp_top_out?.summary_tsv
            def scoring_fade_bot_ch      = params.fade ? FADE.out.summary_bottom : fade_precomp_bot_out?.summary_tsv
            def scoring_fade_site_top_ch = params.fade ? FADE.out.site_tsv_top   : fade_precomp_top_out?.site_tsv
            def scoring_fade_site_bot_ch = params.fade ? FADE.out.site_tsv_bottom: fade_precomp_bot_out?.site_tsv
            def scoring_rer_ch           = (params.rer_tool || params.rer_continuous_file) ? RER_MAIN.out.summary_tsv : null
            def scoring_rer_perms_ch     = (params.rer_tool || params.rer_continuous_file) ? RER_MAIN.out.perms      : null
            def scoring_accum_ch         = accum_results     ? accum_results.results               : null
            def scoring_vep_pai_ch       = params.vep        ? VEP.out.primateai_tsv              : null
            def scoring_vep_cosmic_ch    = params.vep        ? VEP.out.cosmic_tsv                 : null
            // genomic_info comes from params.gene_ensembl_file (resolved inside scoring.nf)

            // CAAS permulation-excess null → genes×N matrices (caas_perms.rds) feeding
            // the CAAS FCS p.perm path. Gated on --caas_permulation_enrichment; needs
            // the full-pool perm-discovery + resample subset emitted by CT.
            def scoring_caas_perms_ch = null
            def scoring_caas_perm_scores_ch = null
            def scoring_caas_pos_pval_ch = null
            def scoring_caas_pos_sample_ch = null
            def caas_perm_out = null
            if (params.caas_permulation_enrichment && ct_results) {
                def caas_universe_ch = (pp_cleaned_bg ?: Channel.empty()).ifEmpty { file('NO_FILE') }
                caas_perm_out = CAAS_PERMULATION(
                    ct_results.caas_perm_discovery,
                    ct_results.caas_resample_subset,
                    ct_results.tree_file,
                    caas_universe_ch
                )
                scoring_caas_perms_ch = caas_perm_out.perms
                scoring_caas_perm_scores_ch = Channel.empty()
                scoring_caas_pos_pval_ch = caas_perm_out.pos_pval      // lean recovery p-value per (gene,position,scheme)
                scoring_caas_pos_sample_ch = caas_perm_out.pos_sample  // lean capped sample for distribution plots
            }

            SCORING(
                scoring_postproc_ch,
                scoring_fade_top_ch,
                scoring_fade_bot_ch,
                scoring_rer_ch,
                scoring_accum_ch,
                scoring_vep_pai_ch,
                scoring_vep_cosmic_ch,
                null,  // genomic_info — resolved from params.gene_ensembl_file in scoring.nf
                scoring_fade_site_top_ch,
                scoring_fade_site_bot_ch,
                pp_cleaned_bg,         // cleaned_background_main.txt — FCS universe
                scoring_rer_perms_ch,  // RER permulation RDS → p.perm in centralized RER FCS
                scoring_caas_perms_ch, // CAAS permulation RDS (asr+caas null) → FCS p.perm + report
                scoring_caas_pos_pval_ch,    // lean recovery p-value per (gene,position,scheme)
                scoring_caas_pos_sample_ch   // lean capped sample for report distribution plots
            )
            ran_any = true

            if (params.enrichment) {
                def fcs_stats_ch = params.scoring ? SCORING.out.fcs_stats : Channel.empty()
                def fcs_stats_rer_ch = params.scoring ? SCORING.out.fcs_stats_rer : Channel.empty()
                def fcs_stats_fade_ch = params.scoring ? SCORING.out.fcs_stats_fade : Channel.empty()
                def fcs_stats_accum_ch = params.scoring ? SCORING.out.fcs_stats_accum : Channel.empty()
                def gene_lists_ch = params.scoring ? SCORING.out.gene_lists : Channel.empty()
                def gene_scores_ch = params.scoring ? SCORING.out.gene_scores : Channel.empty()
                def position_scores_ch = params.scoring ? SCORING.out.position_scores : Channel.empty()
                def position_lists_ch = params.scoring ? SCORING.out.position_lists : Channel.empty()
                // POSENRICH background = caastools background.output (tested positions);
                // the engine restricts it to the cleaned_background genes.
                def posenrich_background_ch = (ct_results && ran_discovery) ? ct_results.background_file : file('NO_FILE')

                // RER's own gene universe + gene lists (significant/accelerating/
                // decelerating) and FADE's per-direction universe + significant
                // gene lists -- feed the unified 13.AMI_analysis.Rmd's FADE/RER
                // sections (see ENRICHMENT workflow). Only populated when the
                // respective tool actually ran this invocation --
                // OR (FADE, mirroring rer_continuous_file's precomputed-input
                // role for RER) when its summary TSV is available from a prior
                // run via --fade_json_dir_top/_bottom without --fade itself
                // this invocation: FADE_GENE_LISTS only needs that summary TSV,
                // so it's re-derived here from fade_precomp_{top,bot}_out the
                // same way scoring_fade_top_ch/scoring_fade_bot_ch already do
                // above, rather than requiring a live --fade run just for AMI.
                def rer_ran  = params.rer_tool || params.rer_continuous_file
                def rer_gene_lists_bg_ch       = rer_ran    ? RER_MAIN.out.gene_lists_bg       : Channel.empty()
                def rer_gene_lists_interest_ch = rer_ran    ? RER_MAIN.out.gene_lists_interest : Channel.empty()

                def fade_gene_lists_bg_top_ch
                def fade_gene_lists_sig_top_ch
                if (params.fade) {
                    fade_gene_lists_bg_top_ch  = FADE.out.gene_lists_bg_top
                    fade_gene_lists_sig_top_ch = FADE.out.gene_lists_sig_top
                } else if (fade_precomp_top_out) {
                    def precomp_top_lists = FADE_GENE_LISTS_PRECOMP_TOP(Channel.value('top'), fade_precomp_top_out.summary_tsv)
                    fade_gene_lists_bg_top_ch  = precomp_top_lists.gene_lists.flatten().filter { it.name == 'background.txt' }.collect()
                    fade_gene_lists_sig_top_ch = precomp_top_lists.gene_lists.flatten().filter { it.name != 'background.txt' }.collect()
                } else {
                    fade_gene_lists_bg_top_ch  = Channel.empty()
                    fade_gene_lists_sig_top_ch = Channel.empty()
                }

                def fade_gene_lists_bg_bottom_ch
                def fade_gene_lists_sig_bottom_ch
                if (params.fade) {
                    fade_gene_lists_bg_bottom_ch  = FADE.out.gene_lists_bg_bottom
                    fade_gene_lists_sig_bottom_ch = FADE.out.gene_lists_sig_bottom
                } else if (fade_precomp_bot_out) {
                    def precomp_bot_lists = FADE_GENE_LISTS_PRECOMP_BOTTOM(Channel.value('bottom'), fade_precomp_bot_out.summary_tsv)
                    fade_gene_lists_bg_bottom_ch  = precomp_bot_lists.gene_lists.flatten().filter { it.name == 'background.txt' }.collect()
                    fade_gene_lists_sig_bottom_ch = precomp_bot_lists.gene_lists.flatten().filter { it.name != 'background.txt' }.collect()
                } else {
                    fade_gene_lists_bg_bottom_ch  = Channel.empty()
                    fade_gene_lists_sig_bottom_ch = Channel.empty()
                }

                ENRICHMENT(
                    fcs_stats_ch,
                    fcs_stats_rer_ch,
                    fcs_stats_fade_ch,
                    fcs_stats_accum_ch,
                    gene_lists_ch,
                    gene_scores_ch,
                    pp_cleaned_bg,
                    scoring_rer_perms_ch,
                    scoring_caas_perms_ch,
                    scoring_caas_perm_scores_ch,
                    scoring_caas_pos_pval_ch,
                    scoring_caas_pos_sample_ch,
                    position_scores_ch,
                    position_lists_ch,
                    posenrich_background_ch,
                    scoring_vep_pai_ch,
                    scoring_vep_cosmic_ch,
                    params.fade ? FADE.out.sites_csv_top    : null,
                    params.fade ? FADE.out.sites_csv_bottom : null,
                    rer_gene_lists_bg_ch,
                    rer_gene_lists_interest_ch,
                    fade_gene_lists_bg_top_ch,
                    fade_gene_lists_bg_bottom_ch,
                    fade_gene_lists_sig_top_ch,
                    fade_gene_lists_sig_bottom_ch
                )
            }
        }

        if (!ran_any) {
            log.info "No tool selected. Use --reporting, --contrast_selection, --ct_tool, --rer_tool, --ct_disambiguation, --ct_postproc, --enrichment, --ct_accumulation, --fade, --scoring, or --rer_tool."
        }
    }
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  THE END: End of the main.nf file.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
