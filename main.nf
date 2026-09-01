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
// Position-level FADE-site CSV (gene,position,max_bf,target_aa) — same process
// FADE() itself calls, re-run standalone against the precomputed *.FADE.json
// dir so posenrich's Position Characterisation FADE-overlap check and
// ENRICHMENT's fade_sites_top/bottom_ch aren't silently null on a
// --fade_json_dir_top/_bottom-only (no live --fade) invocation.
include {FADE_JSON_TO_CSV as FADE_JSON_TO_CSV_PRECOMP_TOP; FADE_JSON_TO_CSV as FADE_JSON_TO_CSV_PRECOMP_BOTTOM} from './subworkflows/FADE/fade_json_to_csv.nf'
include {SELECTION_PREP} from './subworkflows/SELECTION/selection_prep.nf'
include {VEP}            from './workflows/vep.nf'
include {SCORING}        from './workflows/scoring.nf'
include {CAAS_PERMULATION; CAAS_PERMS_PREP} from './subworkflows/CT/caas_permulation.nf'
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
        def fade_precomp_sites_top_ch = null
        def fade_precomp_sites_bot_ch = null

        // FADE (like --ct_tool) pulls in CONTRAST_SELECTION below, which runs
        // REPORTING() itself when --reporting is set. Skip the standalone call
        // here in that case to avoid invoking REPORTING() twice.
        if (params.reporting && !params.contrast_selection && !params.fade) {
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

                def trait_input_for_ct = (contrast_out && contrast_out.trait_dir_out) ? contrast_out.trait_dir_out : contrast_out.trait_file_out
                ct_results = CT(trait_input_for_ct, contrast_out.bootstrap_trait_file_out, contrast_out.tree_file_out)

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

            // Same graceful stop as the CT branch above. This path is reached when
            // CAAStools output is reused but contrast selection still runs to supply
            // the trait file (e.g. recomputing disambiguation over a precomputed
            // discovery). Without it, an under-powered trait emits no traitfile_ok.tab
            // and the run fails further down reporting a missing --caas_config, which
            // says nothing about the actual cause.
            contrast_out.low_contrasts_skip.view { skip_file ->
                exit 0, "Minimum contrast threshold not met for trait '${params.traitname ?: 'unknown'}' (flag: ${skip_file}). Stopping pipeline gracefully."
            }
            ran_any = true
        }
        // FADE needs the foreground/background partition defined by
        // 4.Independent_contrasts.Rmd. When FADE runs standalone (no --ct_tool,
        // no --contrast_selection) trigger CONTRAST_SELECTION here, the same way
        // --ct_tool does — FADE stays independently runnable, the user only
        // passes --fade.
        if (params.fade && !contrast_out) {
            contrast_out = CONTRAST_SELECTION()
            contrast_out.low_contrasts_skip.view { skip_file ->
                exit 0, "Minimum contrast threshold not met for trait '${params.traitname ?: 'unknown'}' (flag: ${skip_file}). Stopping pipeline gracefully."
            }
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

        def toBool = { val ->
            if (val == null) return false
            if (val instanceof Boolean) return val
            if (val instanceof String) return !(val.trim().toLowerCase() in ['false', '0', 'no', 'f', ''])
            return (boolean) val
        }

        def run_ct_disambiguation = toBool(params.ct_disambiguation) && (run_signification || params.signification_from || params.disambiguation_input)
        def run_ct_postproc       = toBool(params.ct_postproc) && (run_ct_disambiguation || params.disambiguation_input)
        def run_ct_accumulation   = toBool(params.ct_accumulation) && (run_ct_postproc || params.accumulation_background_input)
        def run_caas_permulation  = run_ct_disambiguation || (toBool(params.enrichment) && toBool(params.caas_permulation_enrichment))

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

        def scoring_caas_perms_ch = null
        def scoring_caas_perm_scores_ch = null
        def scoring_caas_pos_pval_ch = null
        def scoring_caas_pos_sample_ch = null
        def scoring_caas_pos_quantiles_ch = null
        def caas_perm_out = null

        if (run_ct_disambiguation) {
            // Forward both possible signification metadata artifacts; CT_DISAMBIGUATION
            // will prefer global_meta_caas.tsv when present and otherwise accept
            // the per-run meta_caas.tsv fallback.
            def meta_for_disambiguation = signification_results
                ? signification_results.signification_global_meta.mix(signification_results.signification_meta_caas)
                : null
            // Disambiguation needs the fg/bg trait file(s) that defined the contrasts the
            // CAAS were discovered under. Two suppliers, in order of preference:
            //   • CT ran live            -> ct_results.trait_file (carries traitfiles_ok_dir in multi-hypothesis mode)
            //   • CONTRAST_SELECTION ran -> (contrast_out.trait_dir_out ?: contrast_out.trait_file_out)
            // CONTRAST_SELECTION is deterministic given the same --my_traits and tree,
            // so the pairing it emits here matches the one the precomputed discovery
            // used. Falling through to Channel.empty() lands on --caas_config, which
            // only standalone (non-GUI) runs set. Same resolution as the stats/tree
            // channels built for POSENRICH below.
            def trait_for_disambiguation = ct_results
                ? ct_results.trait_file
                : (contrast_out ? (contrast_out.trait_dir_out ?: contrast_out.trait_file_out) : Channel.empty())
            def tree_for_disambiguation = ct_results
                ? ct_results.tree_file
                : (contrast_out ? contrast_out.tree_file_out : Channel.empty())

            disambiguation_results = CT_DISAMBIGUATION(meta_for_disambiguation, trait_for_disambiguation, tree_for_disambiguation)
            ran_any = true
        }

        if (run_caas_permulation) {


            def perm_disc_ch = null
            def perm_subset_ch = null
            def perm_tree_ch = ct_results ? ct_results.tree_file : (contrast_out ? contrast_out.tree_file_out : (params.tree ? file(params.tree) : Channel.empty()))

            if (ct_results && ct_results.caas_perm_discovery && ct_results.caas_resample_subset) {
                perm_disc_ch = ct_results.caas_perm_discovery
                perm_subset_ch = ct_results.caas_resample_subset
            } else {
                // 1. Check if precomputed permulation discovery outputs already exist from exploratory pass:
                def precomp_disc_files = null
                def precomp_subset_file = null

                def candidate_base_dirs = []
                if (params.resample_from)           candidate_base_dirs.add(file(params.resample_from).parent.parent)
                if (params.disambiguation_input)    candidate_base_dirs.add(file(params.disambiguation_input).parent.parent)
                if (params.discovery_from)         candidate_base_dirs.add(file(params.discovery_from).parent.parent)
                if (params.background_input)        candidate_base_dirs.add(file(params.background_input).parent.parent)
                if (params.signification_from)     candidate_base_dirs.add(file(params.signification_from).parent.parent)
                candidate_base_dirs.add(file(params.outdir))

                for (base_dir in candidate_base_dirs) {
                    if (!base_dir || !base_dir.exists()) continue

                    // 1. Resample file resolution
                    def subset_f = null
                    def resample_candidates = [
                        file("${base_dir}/caas_permulation/resample_perms.tab"),
                        file("${base_dir}/caastools/resample.tab"),
                        file("${base_dir}/caastools/resample"),
                        params.resample_from ? file(params.resample_from) : null
                    ]
                    for (rf in resample_candidates) {
                        if (rf && rf.exists()) { subset_f = rf; break }
                    }

                    // 2. Discovery / Bootstrap permulation files resolution
                    // NOTE: Must be per-cycle *.bootstrap.discovery.output files. Summary bootstrap.tab
                    // lacks the per-cycle 'cycle' column needed by disambiguation_perms_main.py.
                    def disc_candidates = []
                    def disc_dir = file("${base_dir}/caas_permulation/perm_disc")
                    def runtime_boot_dir = file("${base_dir}/runtime/filter/bootstrap")
                    def caas_boot_dir = file("${base_dir}/caastools/bootstrap")

                    if (disc_dir.exists() && file("${disc_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${disc_dir}/*.bootstrap.discovery.output")
                    } else if (runtime_boot_dir.exists() && file("${runtime_boot_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${runtime_boot_dir}/*.bootstrap.discovery.output")
                    } else if (caas_boot_dir.exists() && file("${caas_boot_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${caas_boot_dir}/*.bootstrap.discovery.output")
                    } else if (file("${base_dir}/caas_permulation/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${base_dir}/caas_permulation/*.bootstrap.discovery.output")
                    } else if (file("${base_dir}/caastools/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${base_dir}/caastools/*.bootstrap.discovery.output")
                    }

                    if (disc_candidates && subset_f) {
                        precomp_disc_files = disc_candidates
                        precomp_subset_file = subset_f
                        break
                    }
                }

                if (precomp_disc_files && precomp_subset_file) {
                    log.info "[CAAS_PERMULATION] Reusing precomputed permulation discovery (${precomp_disc_files.size()} file(s)) + resample file (${precomp_subset_file.name}) for ASR re-disambiguation"
                    perm_disc_ch = Channel.fromPath(precomp_disc_files).collect()
                    perm_subset_ch = Channel.value(precomp_subset_file)
                } else {
                    // 2. Precomputed CAAStools run without precomputed perm_disc: check if resample source is available for CAAS_PERMS_PREP
                    def resample_src = null
                    def candidate_resamples = []
                    if (params.resample_from) candidate_resamples.add(params.resample_from)
                    if (params.discovery_from) {
                        def p = file(params.discovery_from).parent
                        candidate_resamples.add("${p}/resample.tab")
                        candidate_resamples.add("${p}/resample")
                    }
                    if (params.background_input) {
                        def p = file(params.background_input).parent
                        candidate_resamples.add("${p}/resample.tab")
                        candidate_resamples.add("${p}/resample")
                    }
                    candidate_resamples.add("${params.outdir}/caastools/resample.tab")
                    candidate_resamples.add("${params.outdir}/caastools/resample")

                    for (cp in candidate_resamples) {
                        if (cp && file(cp).exists()) {
                            resample_src = file(cp)
                            break
                        }
                    }

                    if (resample_src && params.alignment) {
                        def align_dir_f = file(params.alignment)
                        def all_ali_files = align_dir_f.exists() ? (align_dir_f.listFiles()?.findAll { it.isFile() && !it.name.matches('.*\\.txt$|.*\\.tsv$|.*\\.csv$|.*\\.log$|.*\\.map$') } ?: []) : []
                        if (all_ali_files) {
                            if (params.toy_mode) {
                                def n = (params.toy_n ?: 50) as int
                                Collections.shuffle(all_ali_files, new Random((params.seed ?: 1998) as long))
                                all_ali_files = all_ali_files.take(n)
                            }
                            def align_tuple_standalone = Channel.fromList(all_ali_files.collect { f -> tuple(f.baseName, f) })
                            def caas_cfg_standalone = contrast_out ? contrast_out.trait_file_out : (params.caas_config ? file(params.caas_config) : (ct_results ? ct_results.trait_file : Channel.empty()))

                            def perms_prep = CAAS_PERMS_PREP(align_tuple_standalone, caas_cfg_standalone, resample_src)
                            perm_disc_ch = perms_prep.perm_discovery
                            perm_subset_ch = perms_prep.resample_subset
                        }
                    }
                }
            }

            if (perm_disc_ch && perm_subset_ch) {
                def caas_universe_ch = (pp_cleaned_bg ?: Channel.empty()).ifEmpty { file('NO_FILE') }
                // In asr_mode=compute the live CT_DISAMBIGUATION_RUN writes the
                // shared ASR cache that CAAS_PERMS_DISAMBIGUATE reads — gate the
                // replay on it completing (its master_csv is written only after
                // every gene's ASR is cached). No gate when ASR is precomputed
                // (cache already on disk) or disambiguation didn't run live.
                def asr_gate_ch = (params.ct_disambiguation
                                   && (params.ct_disambig_asr_mode ?: 'precomputed') == 'compute'
                                   && disambiguation_results)
                    ? disambiguation_results.master_csv
                    : Channel.value('NO_GATE')
                caas_perm_out = CAAS_PERMULATION(
                    perm_disc_ch,
                    perm_subset_ch,
                    perm_tree_ch,
                    caas_universe_ch,
                    asr_gate_ch
                )
                scoring_caas_perms_ch = caas_perm_out.perms
                scoring_caas_perm_scores_ch = Channel.empty()
                scoring_caas_pos_pval_ch = caas_perm_out.pos_pval      // LOO null_pvalue_boot per (gene,position,scheme)
                scoring_caas_pos_sample_ch = caas_perm_out.pos_sample  // cycle-stratified sample for distribution plots
                scoring_caas_pos_quantiles_ch = caas_perm_out.pos_quantiles  // per (cycle,scheme) distribution shape
                ran_any = true
            }
        }


        if (run_ct_postproc) {
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

        // Fallback resolution for precomputed / standalone runs where --ct_postproc did not run live
        if (!pp_cleaned_bg) {
            def bg_candidate = params.background_input ?: (params.scoring_background_input ?: (params.accumulation_background_input ?: ''))
            if (!bg_candidate && params.scoring_postproc_input) {
                def pfile = file(params.scoring_postproc_input)
                def pdir = pfile ? pfile.parent : null
                if (pdir && file("${pdir}/cleaned_background_main.txt").exists()) {
                    bg_candidate = "${pdir}/cleaned_background_main.txt"
                }
            }
            if (bg_candidate && file(bg_candidate).exists()) {
                pp_cleaned_bg = Channel.value(file(bg_candidate))
            }
        }

        def accum_results = null

        if (run_ct_accumulation) {

            if (!params.ct_postproc && !params.accumulation_background_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone background file (--accumulation_background_input)."
            }
            if (!params.ct_postproc && !params.accumulation_caas_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone CAAS file (--accumulation_caas_input)."
            }

            // Use filtered_discovery.tsv from postproc (gene_filtering stage)
            def acc_caas_ch       = postproc_results ? postproc_results.filtered_discovery : Channel.empty()
            def acc_background_ch = pp_cleaned_bg    ?: Channel.empty()
            def acc_trait_file_ch = ct_results       ? ct_results.trait_file
                : (contrast_out ? (contrast_out.trait_dir_out ?: contrast_out.trait_file_out) : Channel.empty())
            // background.output = the positions CAAStools actually TESTED. This is the
            // accumulation null's eligible pool (intersected with the cleaned-background
            // genes inside the subworkflow). Resolved exactly like POSENRICH's own
            // background: live CT output when discovery ran, else the precomputed param.
            def acc_tested_pos_ch = (ct_results && ran_discovery)
                ? ct_results.background_file
                : (params.posenrich_background_file
                    ? Channel.fromPath(params.posenrich_background_file)
                    : Channel.empty())

            accum_results = CT_ACCUMULATION(acc_caas_ch, acc_background_ch, acc_trait_file_ch,
                                            acc_tested_pos_ch)
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

            // Foreground/background pool for FADE: the PRE-Dunn candidate species
            // from 3.CI-composition.Rmd (candidate_species.tab, traitfile format).
            // FADE tests directional selection on foreground branches and does not
            // need mutually-independent pairs, so it uses the full candidate pool
            // rather than the Dunn-composited canonical set that CAAS uses.
            // contrast_out is always populated here — the block above triggers
            // CONTRAST_SELECTION for standalone --fade.
            def species_source_ch = contrast_out
                ? contrast_out.candidate_species_out
                : Channel.empty()

            // Tree: CT-pruned tree from contrast_selection; otherwise
            // SELECTION_PREP falls back to params.tree via its .ifEmpty {} guard.
            def tree_source_ch = contrast_out
                ? contrast_out.tree_file_out
                : Channel.empty()

            // CT discovery output for toy_mode gene reuse (null when CT didn't run)
            def ct_discovery_source_ch = (ct_results && ran_discovery)
                ? ct_results.discovery_file
                : Channel.empty()

            // gene_set mode sources its directional gene lists from the
            // --fade_postproc_top / --fade_postproc_bottom params (resolved inside
            // SELECTION_PREP); the empty channels here keep the take: signature.
            def sel_pp_top_ch    = Channel.empty()
            def sel_pp_bottom_ch = Channel.empty()

            // Run alignment prep ONCE for FADE.
            // SELECTION_PREP outputs value channels for species/tree files
            // (can be safely consumed by multiple downstream operators) and
            // queue channels for the per-gene filtered FASTAs.
            SELECTION_PREP(
                species_source_ch,
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
            // Position-level FADE-site CSV -- posenrich's Position
            // Characterisation FADE-overlap check needs this (gene,position,
            // max_bf,target_aa), a different file than summary_tsv/site_tsv
            // above (see fade_json_to_csv.nf), so it needs its own precomputed
            // re-run rather than falling out of FADE_REPORT_PRECOMP_*.
            fade_precomp_sites_top_ch = FADE_JSON_TO_CSV_PRECOMP_TOP(Channel.value('top'), precomp_top_jsons).sites_csv
            fade_precomp_sites_bot_ch = FADE_JSON_TO_CSV_PRECOMP_BOTTOM(Channel.value('bottom'), precomp_bottom_jsons).sites_csv
            ran_any = true
        }

        // rer_continuous_file lets RER_MAIN render 5.RERconverge_report.html from a
        // precomputed *.continuous.output RDS even when --rer_tool itself is off
        // this invocation (see rerconverge.nf's own params.rer_continuous_file branch).
        if (params.rer_tool || params.rer_continuous_file) {
            // NOTE: RER_TRAIT requires the original phenotype file (with proper column
            // headers), NOT the caastools traitfile (headerless 3-col format).
            def rer_traitfile_ch = Channel.empty()
            RER_MAIN(
                rer_traitfile_ch,
                Channel.empty(),
                Channel.empty()
            )
            ran_any = true
        }

        println "DEBUG: params.traitname = '${params.traitname}'"

        // CAAS permulation-excess null → genes×N matrices (caas_perms.rds) + the
        // lean position-level LOO null_pvalue_boot. Runs whenever caas_permulation_enrichment
        // is enabled. If live CT ran, consumes ct_results channels; if CT is precomputed
        // (RUN_CAAS=false), resolves precomputed resample + alignment inputs to run
        // CAAS_PERMS_PREP and CAAS_PERMULATION.


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
            // scoring_caas_* are built above, outside this block — see the comment there.

            SCORING(
                scoring_postproc_ch,
                scoring_fade_top_ch,
                scoring_fade_bot_ch,
                scoring_rer_ch,
                scoring_accum_ch,
                null,  // genomic_info — resolved from params.gene_ensembl_file in scoring.nf
                scoring_fade_site_top_ch,
                scoring_fade_site_bot_ch,
                pp_cleaned_bg,         // cleaned_background_main.txt — FCS universe
                scoring_rer_perms_ch,  // RER permulation RDS → p.perm in centralized RER FCS
                scoring_caas_perms_ch, // CAAS permulation RDS (asr+caas null) → FCS p.perm + report
                scoring_caas_pos_pval_ch,    // LOO null_pvalue_boot per (gene,position,scheme)
                scoring_caas_pos_sample_ch,  // cycle-stratified sample for report distribution plots
                scoring_caas_pos_quantiles_ch // per (cycle,scheme) null distribution shape
            )
            ran_any = true

            if (params.enrichment) {
                // SCORING may have REBUILT the CAAS null from --caas_pos_detail_file
                // instead of importing caas_perms.rds. Take the null it actually
                // resolved so ENRICHMENT's FCS p.perm uses the same matrices the
                // scoring report did; re-resolving here would silently fall back to
                // the cached (possibly stale) file.
                def caas_perms_for_enrich = params.scoring ? SCORING.out.caas_perms : scoring_caas_perms_ch
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
                    caas_perms_for_enrich,
                    scoring_caas_perm_scores_ch,
                    scoring_caas_pos_pval_ch,
                    scoring_caas_pos_sample_ch,
                    position_scores_ch,
                    position_lists_ch,
                    posenrich_background_ch,
                    scoring_vep_pai_ch,
                    scoring_vep_cosmic_ch,
                    params.fade ? FADE.out.sites_csv_top    : fade_precomp_sites_top_ch,
                    params.fade ? FADE.out.sites_csv_bottom : fade_precomp_sites_bot_ch,
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
