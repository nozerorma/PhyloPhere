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

version = "1.0.0"

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
 *  NAMED WORKFLOW FOR PIPELINE: This section includes the main CT, ORA and RERConverge workflows.
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
include {ORA} from './workflows/ora.nf'
include {ORA_EXCLUDED} from './workflows/ora_excluded.nf'
include {CT_ACCUMULATION} from './workflows/ct_accumulation.nf'
include {FADE}           from './workflows/fade.nf'
include {MOLERATE}       from './workflows/molerate.nf'
include {PGLS}           from './workflows/pgls.nf'
include {VEP}            from './workflows/vep.nf'
include {SCORING}        from './workflows/scoring.nf'

// Workflow-map helper logic lives in lib/WorkflowMap.groovy (auto-loaded by Nextflow)

// Post-completion: write workflow_map.html again after all publishDir copies are done.
workflow.onComplete {
    try {
        // Resolve to an absolute canonical path so the file is always written
        // to the correct location regardless of JVM working directory at hook time.
        def outdirRaw = params.outdir ? params.outdir.toString() : "${workflow.projectDir}/Out"
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
        def pgls_trait_ch = null
        def pgls_tree_ch = null

        if (params.reporting && !params.contrast_selection) {
            reporting_results = REPORTING()
            pgls_trait_ch = reporting_results.pruned_trait_file
            pgls_tree_ch  = reporting_results.pruned_tree_file
            ran_any = true
        }
        def ct_results
        if (params.ct_tool) {
            if (params.contrast_selection) {
                contrast_out = CONTRAST_SELECTION()
                pgls_trait_ch = contrast_out.pruned_trait_file
                pgls_tree_ch  = contrast_out.pruned_tree_file

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
            pgls_trait_ch = contrast_out.pruned_trait_file
            pgls_tree_ch  = contrast_out.pruned_tree_file
            ran_any = true
        }
        def signification_results = null
        def disambiguation_results = null
        def postproc_results = null

        // Channels populated by CT_POSTPROC when it runs.
        // Used by FADE/MOLERATE/RER gene-set piping without manual path specification.
        def sel_pp_top_ch     = Channel.empty()
        def sel_pp_bottom_ch  = Channel.empty()

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

        // Stable channel references for CT_POSTPROC outputs used by multiple consumers.
        // Populated inside the ct_postproc block when --ct_postproc is enabled.
        def pp_cleaned_bg     = null   // cleaned_background_main (single file, value channel)
        def pp_gene_lists_val = null   // .collect()-ed list of ORA gene list files (value channel)

        if (params.ct_signification) {
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
            if (!params.ct_signification && !params.ct_disambig_caas_metadata) {
                error "CT disambiguation requires signification outputs (--ct_signification) or a standalone metadata file (--ct_disambig_caas_metadata)."
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

        if (params.ora && !params.ct_postproc) {
            // Standalone ORA: no upstream postproc channels — ORA will fall back to
            // --ora_gene_lists_input + --ora_background_input via its .ifEmpty {} guards.
            if (!params.ora_gene_lists_input || !params.ora_background_input) {
                error "Standalone ORA requires --ora_gene_lists_input and --ora_background_input (or enable --ct_postproc for integrated mode)."
            }
            // Pass null for gene_lists so the if(channel) guard in ora.nf detects absence
            // and loads directly from --ora_gene_lists_input (Channel.empty() would silently
            // stall .collect() without ever emitting).
            ORA(Channel.empty(), null)
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
            // Collect ora_gene_lists_files into a value channel so FADE/MOLERATE/RER can each
            // independently filter/flatMap the full list without racing on queue items.
            pp_cleaned_bg     = postproc_results.cleaned_background
            pp_gene_lists_val = postproc_results.ora_gene_lists_files.collect()

            if (params.ora) {
                // ORA on excluded gene lists only: characterises genes removed during
                // post-processing using the original (pre-cleanup) background.
                // Post-proc gene-list ORA is intentionally omitted here; enrichment
                // on significant positions is handled downstream by SCORING (scoring_ora).
                ORA_EXCLUDED(postproc_results.background_ori, postproc_results.excluded_gene_lists_files)
                ran_any = true
            }
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

        if (params.pgls) {
            def pgls_caas_ch = postproc_results ? postproc_results.filtered_discovery : Channel.empty()
            PGLS(pgls_caas_ch, pgls_trait_ch, pgls_tree_ch)
            ran_any = true
        }

        if (params.vep) {
            def vep_caas_ch = postproc_results ? postproc_results.filtered_discovery : Channel.empty()
            VEP(vep_caas_ch)
            ran_any = true
        }

        // Populate postproc gene-list channels for FADE/MOLERATE/RER gene-set piping.
        // pp_gene_lists_val is a value channel holding the collected List<File>.
        // Restrict to global scenario directional lists only:
        //   all_top.txt and all_bottom.txt
        // (postproc-only selection policy requested for FADE/MOLERATE).
        if (pp_gene_lists_val) {
            sel_pp_top_ch    = pp_gene_lists_val
                .flatMap { files -> files.findAll { f ->
                    f.name == 'all_top.txt'
                } }
            sel_pp_bottom_ch = pp_gene_lists_val
                .flatMap { files -> files.findAll { f ->
                    f.name == 'all_bottom.txt'
                } }
        }

        // Merge per-phenotype files into single files before passing to selection
        // workflows. When multiple phenotypes are run together each sel_* channel
        // emits N items (one per phenotype). COLLECT_GENE_SETS inside FADE/MOLERATE/RER
        // is invoked once per synchronised 2-tuple of inputs, so without this merge it
        // runs N times and every gene appears N times in the resulting ali channel —
        // causing file-name staging collisions in the report process. collectFile()
        // collapses N items into one merged file so COLLECT_GENE_SETS runs exactly once.
        // When sel_* channels are empty (standalone runs) collectFile() emits nothing
        // and the .ifEmpty {} fallback paths inside each workflow remain intact.
        sel_pp_top_ch     = sel_pp_top_ch    .collectFile(name: 'merged_pp_top.txt')
        sel_pp_bottom_ch  = sel_pp_bottom_ch .collectFile(name: 'merged_pp_bottom.txt')

        // Split stats channel so FADE and MOLERATE can each receive the single
        // emission without competing on the same queue channel.
        def stats_source_ch = contrast_out
            ? contrast_out.stats_file_out
            : (reporting_results ? reporting_results.stats_file : Channel.empty())
        def split_stats_file = stats_source_ch
            ? stats_source_ch.multiMap { f -> fade: f; molerate: f }
            : Channel.empty().multiMap { f -> fade: f; molerate: f }
        // Tree: FADE and MOLERATE always use params.tree (the complete species tree)
        // via their .ifEmpty {} fallback.  The CT-pruned tree is pruned to
        // valid-contrast-pair species only, which is a strict subset of the
        // extreme-labeled species in stats_file — passing it would cause
        // ANNOTATE_TREE_FG / MOLERATE_RUN to silently drop FG branches.

        if (params.fade) {
            FADE(split_stats_file.fade, Channel.empty(),
                 sel_pp_top_ch,  sel_pp_bottom_ch)
            ran_any = true
        }

        if (params.molerate) {
            MOLERATE(split_stats_file.molerate, Channel.empty(),
                     sel_pp_top_ch,  sel_pp_bottom_ch)
            ran_any = true
        }

        if (params.rer_tool) {
            // RER_MAIN is called last so that sel_pp_* channels are fully
            // populated by CT_POSTPROC when running together.
            // NOTE: RER_TRAIT requires the original phenotype file (with proper column
            // headers), NOT the caastools traitfile (headerless 3-col format).
            // Always pass Channel.empty() so RER_MAIN falls back to --my_traits.
            def rer_traitfile_ch = Channel.empty()
            RER_MAIN(
                rer_traitfile_ch,
                sel_pp_top_ch,  sel_pp_bottom_ch
            )
            ran_any = true
        }

        if (params.scoring) {
            if (!params.ct_postproc && !params.scoring_postproc_input) {
                error "SCORING requires CT post-processing output (--ct_postproc) or --scoring_postproc_input."
            }

            // Wire upstream outputs into SCORING. Pass null (not Channel.empty())
            // when a module didn't run so if(channel) guards detect absence correctly.
            def scoring_postproc_ch  = postproc_results ? postproc_results.filtered_discovery : null
            def scoring_pgls_ch      = params.pgls      ? PGLS.out.pgls_tsv                  : null
            def scoring_fade_top_ch  = params.fade       ? FADE.out.summary_top               : null
            def scoring_fade_bot_ch  = params.fade       ? FADE.out.summary_bottom            : null
            def scoring_rer_ch       = params.rer_tool   ? RER_MAIN.out.summary_tsv           : null
            def scoring_accum_ch     = accum_results     ? accum_results.results               : null
            def scoring_bg_ch        = pp_cleaned_bg     ?: null
            def scoring_vep_tv_ch    = params.vep        ? VEP.out.transvar_tsv               : null
            def scoring_vep_pai_ch   = params.vep        ? VEP.out.primateai_tsv              : null
            // genomic_info comes from params.gene_ensembl_file (resolved inside scoring.nf)

            SCORING(
                scoring_postproc_ch,
                scoring_pgls_ch,
                scoring_fade_top_ch,
                scoring_fade_bot_ch,
                scoring_rer_ch,
                scoring_accum_ch,
                scoring_bg_ch,
                scoring_vep_tv_ch,
                scoring_vep_pai_ch,
                null  // genomic_info — resolved from params.gene_ensembl_file in scoring.nf
            )
            ran_any = true
        }

        if (!ran_any) {
            log.info "No tool selected. Use --reporting, --contrast_selection, --ct_tool, --rer_tool, --ct_signification, --ct_disambiguation, --ct_postproc, --ora, --ct_accumulation, --fade, --molerate, --pgls, --scoring, or --rer_tool."
        }
    }
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  THE END: End of the main.nf file.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
