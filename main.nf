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

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  NAMED WORKFLOW FOR PIPELINE: This section includes the main workflows.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include {HELP} from './workflows/help.nf'
include {CT} from './workflows/ct.nf'
include {listAlignmentFiles; sampleAlignmentFiles} from './subworkflows/CT/ct_alignment_files.nf'
include {RER_MAIN} from './workflows/rerconverge.nf'
include {REPORTING} from './workflows/reporting.nf'
include {CONTRAST_SELECTION} from './workflows/contrast_selection.nf'
include {CT_META_CAAS} from './workflows/ct_meta_caas.nf'
include {CT_POSTPROC} from './workflows/ct_postproc.nf'
include {CT_OBSERVED} from './workflows/ct_observed.nf'
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
include {VEP}                       from './workflows/vep.nf'
include {SCORING}        from './workflows/scoring.nf'
include {CAAS_SIGNIFICANCE_REPORT} from './subworkflows/CT_META_CAAS/ctpp_meta_caas.nf'
include {CAAS_CORE; CAAS_CORE_OBSERVED; CAAS_CORE_MERGE; CAAS_PERMS_PREP} from './subworkflows/CT/caas_permulation.nf'
include {CAAS_EVIDENCE}   from './subworkflows/CT_DISAMBIGUATION/ct_evidence.nf'
include {ENRICHMENT}      from './workflows/enrichment.nf'

// Coerce a param that may arrive as Boolean, String ("false", "0", ...) or null into a Boolean.
def toBool(val) {
    if (val == null) return false
    if (val instanceof Boolean) return val
    if (val instanceof String) return !(val.trim().toLowerCase() in ['false', '0', 'no', 'f', ''])
    return val as boolean
}

// Foreground species list published next to a precomputed FADE json dir
// (<selection>/species_sets/<filename>); falls back to the NO_FG_LIST sentinel.
def resolve_fg_species(json_dir, filename) {
    if (!json_dir) return file('NO_FG_LIST')
    def selection_dir = file(json_dir).parent?.parent?.parent
    def candidate = selection_dir ? selection_dir.resolve("species_sets/${filename}") : null
    return (candidate && candidate.exists()) ? candidate : file('NO_FG_LIST')
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RUN PHYLOPHERE ANALYSIS: This section initiates the main Phylophere workflow.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow {

    def version = "2.0.0"

    // Display input parameters
    log.info """

PHYLOPHERE - NF PIPELINE  ~  version ${version}
=============================================

PHYLOPHERE: A Nextflow pipeline including a complete set
of phylogenetic comparative tools and analyses for Phenome-Genome studies

Author:         Miguel Ramon (miguel.ramon@upf.edu)


"""

    // Post-completion: write workflow_map.html again after all publishDir copies are done.
    // WorkflowMap lives in lib/WorkflowMap.groovy (compiled and loaded by Nextflow).
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

    // Check if --help is provided
    if (params.help) {
        HELP ()
    } else {
        // tax_id / gene_ensembl_file auto-generation (bin/resolve_core_inputs.py)
        // happens BEFORE this pipeline is invoked, not here: Nextflow enforces
        // single-assignment on params keys, so a params.tax_id = ... here would
        // be silently ignored once conf/common.config's own params.tax_id = ""
        // default has already run ("`params.tax_id` is defined
        // multiple times -- Assignments following the first are ignored",
        // confirmed empirically). See bin/resolve_core_inputs.py's docstring;
        // the GUI's generated run scripts call it automatically.
        if (!params.tax_id) {
            log.warn "params.tax_id is empty. To auto-generate it, run bin/resolve_core_inputs.py before this pipeline (see its docstring) — it cannot be generated from inside main.nf."
        }
        if (!params.gene_ensembl_file) {
            log.warn "params.gene_ensembl_file is empty. To auto-generate it, run bin/resolve_core_inputs.py before this pipeline (see its docstring) — it cannot be generated from inside main.nf."
        }

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
                ct_results = CT(trait_input_for_ct, contrast_out.permulation_trait_file_out, contrast_out.tree_file_out)

            } else {
                def trait_file_in = null
                def permulation_trait_file_in = null
                def tree_file_in = null
                ct_results = CT (trait_file_in, permulation_trait_file_in, tree_file_in)
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
        // passes --fade. Skipped entirely when --fade_species_file supplies the
        // fg/bg partition directly (candidate_species.tab format), since that
        // bypasses PSS/Dunn contrast selection for FADE.
        if (params.fade && !contrast_out && !params.fade_species_file) {
            contrast_out = CONTRAST_SELECTION()
            contrast_out.low_contrasts_skip.view { skip_file ->
                exit 0, "Minimum contrast threshold not met for trait '${params.traitname ?: 'unknown'}' (flag: ${skip_file}). Stopping pipeline gracefully."
            }
            ran_any = true
        }
        def meta_caas_results = null
        def observed_results = null   // master_csv and results_dir (ct_disambiguation/) of the observed labeling
        def observed_meta = null      // CT_OBSERVED outputs, when a discovery.tab was scored there
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

        // The observed labeling is the b_0 slice of the permulation core: when the core replays the alignments
        // (ct_tool 'discovery', or the standalone replay), CAAS_CORE_OBSERVED writes discovery.tab, the background
        // files, the meta_caas tables and the master CSV from it. A discovery.tab given with --discovery_from is
        // scored by CT_OBSERVED instead, with the same code.
        def run_ct_disambiguation = toBool(params.ct_disambiguation) && (ran_discovery || params.discovery_from || params.disambiguation_input)
        def run_ct_postproc       = toBool(params.ct_postproc) && (run_ct_disambiguation || params.disambiguation_input)
        def run_ct_accumulation   = toBool(params.ct_accumulation) && (run_ct_postproc || params.accumulation_background_input)
        def run_caas_permulation  = run_ct_disambiguation || ran_discovery || (toBool(params.enrichment) && toBool(params.caas_permulation_enrichment))

        // Stable channel references for CT_POSTPROC outputs used by multiple consumers.
        // Populated inside the ct_postproc block when --ct_postproc is enabled.
        def pp_cleaned_bg     = null   // cleaned_background_main (single file, value channel)

        def scoring_caas_perms_ch = null
        def scoring_caas_pos_cycle_caas_ch = null   // perm_pos_cycle_caas.tsv.gz — p.emp numerator/denominator
        def scoring_caas_pos_sample_ch = null
        def scoring_caas_pos_quantiles_ch = null
        def scoring_caas_pos_detail_ch = null        // sharded perm_pos_detail dir — CT_ACCUMULATION permulation null (Tier 3E)
        def scoring_caas_gene_cycle_scores_ch = null // gene_cycle_scores.tsv — CT_ACCUMULATION permulation null (Tier 3E)
        def core = null               // the permulation core: shard directories and b_0 slices of the batches
        def core_observed = null      // the observed contract files, written from the b_0 slices
        def core_replays = false      // true when the core replays the alignments (it then has a b_0 slice)
        def evidence_inputs = null    // [discovery, design, tree] of the observed labeling, for CAAS_EVIDENCE
        def caas_gene_lengths_ch = Channel.value(
            params.gene_ensembl_file ? file(params.gene_ensembl_file) : file('NO_FILE'))
        if (run_caas_permulation) {

            def perm_align_ch = Channel.empty()   // alignments to replay
            def perm_reuse_ch = Channel.empty()   // perm-discovery exports that already exist
            def perm_cfg_ch = Channel.value(file('NO_CONFIG'))
            def perm_subset_ch = null
            def perm_fop_pairs_ch = Channel.value(file('NO_FOP_PAIRS'))
            // Must stay a channel: CAAS_CORE turns it into a value channel with collect().
            // ct_results.tree_file / contrast_out.tree_file_out are workflow emits (already
            // channels); the bare params.tree path needs wrapping or the call throws
            // MissingMethodException on sun.nio.fs.UnixPath.
            def perm_tree_ch = ct_results ? ct_results.tree_file : (contrast_out ? contrast_out.tree_file_out : (params.tree ? Channel.value(file(params.tree)) : Channel.empty()))

            if (ct_results && ct_results.caas_align_tuple && ct_results.caas_resample_subset) {
                perm_align_ch = ct_results.caas_align_tuple
                core_replays = true
                perm_subset_ch = ct_results.caas_resample_subset
                perm_cfg_ch = ct_results.trait_file
                if (ct_results.caas_fop_pairs) perm_fop_pairs_ch = ct_results.caas_fop_pairs
            } else {
                // 1. Check if precomputed permulation discovery outputs already exist from exploratory pass:
                def precomp_disc_files = null
                def precomp_subset_file = null

                def candidate_base_dirs = []
                if (params.resample_from)           candidate_base_dirs.add(file(params.resample_from).parent.parent)
                if (params.disambiguation_input)    candidate_base_dirs.add(file(params.disambiguation_input).parent.parent)
                if (params.discovery_from)         candidate_base_dirs.add(file(params.discovery_from).parent.parent)
                if (params.background_input)        candidate_base_dirs.add(file(params.background_input).parent.parent)
                if (params.meta_caas_from)          candidate_base_dirs.add(file(params.meta_caas_from).parent.parent)
                candidate_base_dirs.add(file(params.outdir))

                // First candidate dir holding both a resample file and per-cycle discovery files wins.
                def precomp_found = candidate_base_dirs.findResult { base_dir ->
                    if (!base_dir || !base_dir.exists()) return null

                    // 1. Resample file resolution
                    def resample_candidates = [
                        file("${base_dir}/caas_permulation/resample_perms.tab"),
                        file("${base_dir}/caastools/resample.tab"),
                        file("${base_dir}/caastools/resample"),
                        params.resample_from ? file(params.resample_from) : null
                    ]
                    def subset_f = resample_candidates.find { rf -> rf && rf.exists() }

                    // 2. Discovery / perm-replay permulation files resolution
                    // NOTE: Must be per-cycle *.perm_replay.discovery.output files (legacy runs:
                    // *.bootstrap.discovery.output). Summary bootstrap.tab lacks the per-cycle
                    // 'cycle' column needed by disambiguation_perms_main.py.
                    def disc_candidates = []
                    def disc_dir = file("${base_dir}/caas_permulation/perm_disc")
                    def runtime_boot_dir = file("${base_dir}/runtime/filter/bootstrap")
                    def caas_boot_dir = file("${base_dir}/caastools/bootstrap")

                    if (disc_dir.exists() && file("${disc_dir}/*.perm_replay.discovery.output")) {
                        disc_candidates = file("${disc_dir}/*.perm_replay.discovery.output")
                    } else if (disc_dir.exists() && file("${disc_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${disc_dir}/*.bootstrap.discovery.output")
                    } else if (runtime_boot_dir.exists() && file("${runtime_boot_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${runtime_boot_dir}/*.bootstrap.discovery.output")
                    } else if (caas_boot_dir.exists() && file("${caas_boot_dir}/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${caas_boot_dir}/*.bootstrap.discovery.output")
                    } else if (file("${base_dir}/caas_permulation/*.perm_replay.discovery.output")) {
                        disc_candidates = file("${base_dir}/caas_permulation/*.perm_replay.discovery.output")
                    } else if (file("${base_dir}/caas_permulation/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${base_dir}/caas_permulation/*.bootstrap.discovery.output")
                    } else if (file("${base_dir}/caastools/*.bootstrap.discovery.output")) {
                        disc_candidates = file("${base_dir}/caastools/*.bootstrap.discovery.output")
                    }

                    return (disc_candidates && subset_f) ? [disc: disc_candidates, subset: subset_f] : null
                }
                precomp_disc_files = precomp_found?.disc
                precomp_subset_file = precomp_found?.subset

                if (precomp_disc_files && precomp_subset_file) {
                    // The exports were made from the labelings file next to them; b_0 is its first row when the
                    // run replayed the real labeling. Exports without it carry no b_0 slice.
                    def first_labeling = precomp_subset_file.withReader { r -> r.readLine() } ?: ''
                    if (!(first_labeling ==~ /^b_0(~[^\t]*)?\t.*/)) {
                        error "[CAAS_CORE] The perm-discovery exports to reuse were made without the real labeling (b_0): ${precomp_subset_file} has no b_0 rows. Rerun the null to replay it."
                    }
                    log.info "[CAAS_CORE] Reusing precomputed permulation discovery (${precomp_disc_files.size()} file(s)) + resample file (${precomp_subset_file.name}) for ASR re-disambiguation"
                    perm_reuse_ch = Channel.fromPath(precomp_disc_files).collect()
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

                    def resample_hit = candidate_resamples.find { cp -> cp && file(cp).exists() }
                    if (resample_hit) resample_src = file(resample_hit)

                    if (resample_src && params.alignment) {
                        def all_ali_files = listAlignmentFiles(params.alignment)
                        if (all_ali_files) {
                            if (params.toy_mode) {
                                def n = (params.toy_n ?: 50) as int
                                all_ali_files = sampleAlignmentFiles(all_ali_files, n, params.seed ?: 1998)
                            }
                            def align_tuple_standalone = Channel.fromList(all_ali_files.collect { f -> tuple(f.baseName, f) })
                            def caas_cfg_standalone = contrast_out ? contrast_out.trait_file_out : (params.caas_config ? file(params.caas_config) : (ct_results ? ct_results.trait_file : Channel.empty()))

                            def perms_prep = CAAS_PERMS_PREP(caas_cfg_standalone, resample_src)
                            perm_align_ch = align_tuple_standalone
                            core_replays = true
                            perm_cfg_ch = (caas_cfg_standalone instanceof java.nio.file.Path) ? Channel.value(caas_cfg_standalone) : caas_cfg_standalone
                            perm_subset_ch = perms_prep.resample_subset
                            perm_fop_pairs_ch = perms_prep.fop_pairs
                        }
                    }
                }
            }

            if (perm_subset_ch) {
                core = CAAS_CORE(
                    perm_align_ch,
                    perm_reuse_ch,
                    perm_cfg_ch,
                    perm_subset_ch,
                    perm_tree_ch,
                    perm_fop_pairs_ch,
                    caas_gene_lengths_ch
                )
                // The b_0 slice of a replay is the observed labeling. A discovery.tab given with --discovery_from
                // is scored by CT_OBSERVED instead, so the two never write the same files.
                if (core_replays && !params.discovery_from && (ran_discovery || run_ct_disambiguation)) {
                    core_observed = CAAS_CORE_OBSERVED(core.b0_observed.collect(), perm_cfg_ch.collect().map { items -> items[0] })
                    evidence_inputs = [discovery: core_observed.discovery, design: perm_cfg_ch, tree: perm_tree_ch]
                }
                ran_any = true
            }
        }

        // The pattern-annotation report reads the observed discovery: the core's, or the one given with --discovery_from.
        def run_meta_caas = (core_observed != null) || params.discovery_from
        if (run_meta_caas) {
            // Only pass the core's channels when it produced the observed files.
            // Pass null (not Channel.empty()) when absent so the if(channel) guard
            // inside CT_META_CAAS correctly detects absence and falls back to params.
            def discovery_ch        = core_observed ? core_observed.discovery        : null
            def background_genes_ch = core_observed ? core_observed.background_genes : null

            meta_caas_results = CT_META_CAAS(discovery_ch, background_genes_ch)
            ran_any = true
        }

        if (run_ct_disambiguation) {
            if (core_observed) {
                observed_results = [master_csv: core_observed.master_csv, results_dir: core_observed.results_dir]
            } else if (params.discovery_from) {
                // A discovery.tab that already exists is scored with the code of the core's b_0 slice.
                // Disambiguation needs the fg/bg trait file(s) that defined the contrasts the
                // CAAS were discovered under. Two suppliers, in order of preference:
                //   • CT ran live            -> ct_results.trait_file (carries traitfiles_ok_dir in multi-hypothesis mode)
                //   • CONTRAST_SELECTION ran -> (contrast_out.trait_dir_out ?: contrast_out.trait_file_out)
                // CONTRAST_SELECTION is deterministic given the same --my_traits and tree,
                // so the pairing it emits here matches the one the precomputed discovery
                // used. Falling through to Channel.empty() lands on --caas_config, which
                // only standalone (non-GUI) runs set.
                def trait_for_observed = ct_results
                    ? ct_results.trait_file
                    : (contrast_out ? (contrast_out.trait_dir_out ?: contrast_out.trait_file_out) : Channel.empty())
                def tree_for_observed = ct_results
                    ? ct_results.tree_file
                    : (contrast_out ? contrast_out.tree_file_out : Channel.empty())

                // FOP pair weights from the same trait supplier: in multi-hypothesis
                // mode the trait channel is a directory (traitfiles_ok_dir /
                // Traitfiles) holding contrast_hypotheses_pairs.tsv.
                def hyp_pairs_for_observed = (ct_results || contrast_out)
                    ? trait_for_observed.map { d ->
                          def f = file("${d}/contrast_hypotheses_pairs.tsv")
                          (file(d).isDirectory() && f.exists()) ? f : file('NO_HYP_PAIRS')
                      }
                    : null

                def discovery_file_obj = file(params.discovery_from)
                assert discovery_file_obj.exists() : "Error: discovery_from file not found: ${params.discovery_from}"
                def observed_run = CT_OBSERVED(Channel.value(discovery_file_obj), trait_for_observed, tree_for_observed, hyp_pairs_for_observed)
                observed_results = [master_csv: observed_run.master_csv, results_dir: observed_run.results_dir]
                observed_meta = observed_run
                evidence_inputs = [discovery: Channel.value(discovery_file_obj), design: observed_run.design, tree: observed_run.tree]
                ran_any = true
            }
        }
        // CT_POSTPROC is resolved before the permulation null so that the null receives the real
        // cleaned background as its gene universe (dataflow, not call order, decides when tasks run).
        if (run_ct_postproc) {
            // Post-processing is downstream from disambiguation; consume disambiguation master CSV when available
            // Pass null (not Channel.empty()) when there is no upstream result so that the
            // if(channel) guard inside CT_POSTPROC correctly detects absence and falls back
            // to --disambiguation_input / --background_input params (same pattern as CT_META_CAAS).
            def disambiguation_ch = observed_results ? observed_results.master_csv : null
            // Only wire the background genes when discovery actually ran; otherwise pass null so
            // CT_POSTPROC falls back to the --background_input param.
            def background_genes_ch = core_observed ? core_observed.background_genes : null
            // Pass full ct_disambiguation/ directory for ASR robustness diagnostics (null = standalone mode)
            def disambiguation_dir_ch = observed_results ? observed_results.results_dir : null
            // Contrast design (top / bottom species per hypothesis) for the species tally of the position table
            def postproc_hyp_pairs_ch = contrast_out
                ? (contrast_out.trait_dir_out ?: Channel.empty())
                      .map { d ->
                          def f = d ? file("${d}/contrast_hypotheses_pairs.tsv") : null
                          (f && f.exists()) ? f : file('NO_HYP_PAIRS')
                      }
                : null
            postproc_results = CT_POSTPROC(disambiguation_ch, background_genes_ch, disambiguation_dir_ch, postproc_hyp_pairs_ch)
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

        if (core) {
            def caas_universe_ch = (pp_cleaned_bg ?: Channel.empty()).ifEmpty { file('NO_FILE') }
            def caas_perm_out = CAAS_CORE_MERGE(core.pos_detail.collect(), caas_universe_ch, caas_gene_lengths_ch.collect().map { items -> items[0] }, core.labelings)
            scoring_caas_perms_ch = caas_perm_out.perms
            scoring_caas_pos_cycle_caas_ch = caas_perm_out.pos_cycle_caas.ifEmpty(file('NO_CAAS_POS_CYCLE_CAAS'))  // per (gene,position,side,cycle) caas_score -> p.emp
            scoring_caas_pos_sample_ch = caas_perm_out.pos_sample.ifEmpty(file('NO_CAAS_POS_SAMPLE'))  // cycle-stratified sample for distribution plots
            scoring_caas_pos_quantiles_ch = caas_perm_out.pos_quantiles.ifEmpty(file('NO_CAAS_POS_QUANTILES'))  // per (cycle,scheme) distribution shape
            scoring_caas_pos_detail_ch = caas_perm_out.pos_detail                   // sharded perm_pos_detail dir
            scoring_caas_gene_cycle_scores_ch = caas_perm_out.gene_cycle_scores     // genes x cycles raw scores
            ran_any = true
        }


        def accum_results = null

        if (run_ct_accumulation) {

            if (!params.ct_postproc && !params.accumulation_background_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone background file (--accumulation_background_input)."
            }
            if (!params.ct_postproc && !params.accumulation_caas_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone CAAS file (--accumulation_caas_input)."
            }

            // The permulation null of this run has no cycle when only b_0 is replayed (--caas_full_perms 0). Written without
            // `?:`: Groovy reads a numeric 0 as false.
            if (core && params.accumulation_randomization_type == 'permulation' && params.caas_full_perms != null && (params.caas_full_perms as int) == 0) {
                error "CT_ACCUMULATION with --accumulation_randomization_type permulation needs the permuted null, which is empty with --caas_full_perms 0: use naive or cons_decile, or raise --caas_full_perms."
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
            def acc_tested_pos_ch = core_observed
                ? core_observed.background
                : (params.posenrich_background_file
                    ? Channel.fromPath(params.posenrich_background_file)
                    : Channel.empty())

            // sharded perm_pos_detail dir + its companion gene_cycle_scores.tsv —
            // only consumed when accumulation_randomization_type == 'permulation'
            // (Tier 3E); NO_FILE sentinel otherwise, with the standalone
            // --caas_pos_detail_file / --caas_gene_cycle_scores_file params as a
            // fallback for reruns without a live permulation core. The gene_cycle_scores
            // file is what gives the permulation null its exact cycle count instead of
            // inferring it from perm_pos_detail alone (see randomize.py).
            def acc_pos_detail_ch = (scoring_caas_pos_detail_ch ?: Channel.empty())
                .ifEmpty { file(params.caas_pos_detail_file ?: 'NO_FILE') }
            def acc_gene_cycle_scores_ch = (scoring_caas_gene_cycle_scores_ch ?: Channel.empty())
                .ifEmpty { file(params.caas_gene_cycle_scores_file ?: 'NO_FILE') }

            accum_results = CT_ACCUMULATION(acc_caas_ch, acc_background_ch, acc_trait_file_ch,
                                            acc_tested_pos_ch, acc_pos_detail_ch, acc_gene_cycle_scores_ch)
            ran_any = true

        }

        // VEP is invoked further down after SCORING: it consumes position_scores.tsv
        // directly from SCORING (SCORING -> VEP -> ENRICHMENT).

        if (params.fade) {
            // Resolve upstream channel sources for SELECTION_PREP.
            // These channels are now consumed by a single SELECTION_PREP call
            // (instead of being split separately for FADE).

            // Foreground/background pool for FADE: either a user-supplied
            // --fade_species_file (same candidate_species.tab format: species
            // <TAB> contrast_group [<TAB> pair, ignored] — contrast_group==1
            // -> top/fg-eligible, ==0 -> bottom/fg-eligible — parsed as-is by
            // EXTRACT_EXTREME_SPECIES/extract_extreme_species.py, no new
            // parser needed), or the PRE-Dunn candidate species from
            // 3.CI-composition.Rmd. FADE tests directional selection on
            // foreground branches and does not need mutually-independent
            // pairs, so it uses the full candidate pool rather than the
            // Dunn-composited canonical set that CAAS uses. contrast_out is
            // null when --fade_species_file bypassed CONTRAST_SELECTION above.
            def species_source_ch = params.fade_species_file
                ? Channel.fromPath(params.fade_species_file, checkIfExists: true)
                : (contrast_out ? contrast_out.candidate_species_out : Channel.empty())

            // Tree: CT-pruned tree from contrast_selection; otherwise
            // SELECTION_PREP falls back to params.tree via its .ifEmpty {} guard
            // (required when --fade_species_file is used standalone, since no
            // CT-pruned tree exists in that path).
            def tree_source_ch = contrast_out
                ? contrast_out.tree_file_out
                : Channel.empty()

            // CT discovery output for toy_mode gene reuse (null when CT didn't run)
            def ct_discovery_source_ch = core_observed
                ? core_observed.discovery
                : Channel.empty()

            // gene_set mode sources its directional gene lists from the
            // the fade precomputed-summary params (resolved inside
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
            // The GUI's Precomputed Run tab points fade_json_dir_{top,bottom} at
            // <outdir>/selection/fade/<direction>/json (see run_single.sh.j2's
            // PRECOMP_OUTDIR block). EXTRACT_EXTREME_SPECIES published that same
            // prior run's foreground list as a sibling under the same selection/
            // root: <outdir>/selection/species_sets/<direction>_species.txt.
            // Without this, 6.FADE_report.Rmd always got the NO_FG_LIST sentinel
            // here (unconditionally, unlike the live path a few lines above,
            // which wires SELECTION_PREP's channel straight through) and printed
            // "Foreground species list not provided to this report" even though
            // the source run had one. Best-effort: falls back to the sentinel
            // when the derived path doesn't resolve (e.g. --fade_json_dir_top
            // pointed somewhere outside this layout).
            if (params.fade_json_dir_top) {
                def precomp_top_jsons = Channel.fromPath("${params.fade_json_dir_top}/*.FADE.json").collect().ifEmpty([])
                def fg_top_precomp    = resolve_fg_species(params.fade_json_dir_top, 'top_species.txt')
                fade_precomp_top_out      = FADE_REPORT_PRECOMP_TOP(Channel.value('top'), precomp_top_jsons, fg_top_precomp)
                fade_precomp_sites_top_ch = FADE_JSON_TO_CSV_PRECOMP_TOP(Channel.value('top'), precomp_top_jsons).sites_csv
            }

            if (params.fade_json_dir_bottom) {
                def precomp_bottom_jsons = Channel.fromPath("${params.fade_json_dir_bottom}/*.FADE.json").collect().ifEmpty([])
                def fg_bottom_precomp    = resolve_fg_species(params.fade_json_dir_bottom, 'bottom_species.txt')
                fade_precomp_bot_out      = FADE_REPORT_PRECOMP_BOTTOM(Channel.value('bottom'), precomp_bottom_jsons, fg_bottom_precomp)
                fade_precomp_sites_bot_ch = FADE_JSON_TO_CSV_PRECOMP_BOTTOM(Channel.value('bottom'), precomp_bottom_jsons).sites_csv
            }
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
        // per-cycle position-level null (perm_pos_cycle_caas.tsv.gz, feeds p.emp).
        // Runs whenever caas_permulation_enrichment is enabled. If live CT ran,
        // consumes ct_results channels; if CT is precomputed
        // (RUN_CAAS=false), resolves precomputed resample + alignment inputs to run
        // CAAS_PERMS_PREP and CAAS_CORE.


        def evidence_top_n = (params.caas_evidence_top_n ?: 0) as int
        if (evidence_top_n > 0 && !params.scoring) {
            error "caas_evidence_top_n > 0 explains the best positions of position_scores.tsv: it needs --scoring."
        }
        if (params.scoring) {
            if (!params.ct_postproc && !params.scoring_postproc_input) {
                error "SCORING requires CT post-processing output (--ct_postproc) or --scoring_postproc_input."
            }

            // Wire upstream outputs into SCORING. Pass null (not Channel.empty())
            // when a module didn't run so if(channel) guards detect absence correctly.
            def scoring_postproc_ch      = postproc_results
                ? postproc_results.filtered_discovery.collect().map { it[0] }
                : (params.scoring_postproc_input ? Channel.value(file(params.scoring_postproc_input)) : null)
            // Precomputed FADE JSONs (--fade_json_dir_top/_bottom, no --fade this run)
            // feed the exact same summary_tsv/site_tsv SCORING needs, via
            // fade_precomp_{top,bot}_out captured above — so a --scoring run against a
            // prior --fade run's JSON output no longer needs separately pre-built
            // --scoring_fade_summary_*/--scoring_fade_site_* TSVs.
            def scoring_fade_top_ch      = params.fade ? FADE.out.summary_top    : fade_precomp_top_out?.summary_tsv
            def scoring_fade_bot_ch      = params.fade ? FADE.out.summary_bottom : fade_precomp_bot_out?.summary_tsv
            // Position-level FADE evidence for scoring: the per-significant-site
            // table (gene,position,max_bf,target_aa) from FADE_JSON_TO_CSV /
            // parse_fade_json_sites.R — NOT the report's old site_tsv emit,
            // which 6.FADE_report.Rmd no longer produces. scoring_compute.R's
            // .load_fade_sites() is written for exactly this comma-delimited
            // schema, and ENRICHMENT already consumes the same channel.
            def scoring_fade_site_top_ch = params.fade ? FADE.out.sites_csv_top    : fade_precomp_sites_top_ch
            def scoring_fade_site_bot_ch = params.fade ? FADE.out.sites_csv_bottom : fade_precomp_sites_bot_ch
            def scoring_rer_ch           = (params.rer_tool || params.rer_continuous_file) ? RER_MAIN.out.summary_tsv : null
            def scoring_rer_perms_ch     = (params.rer_tool || params.rer_continuous_file) ? RER_MAIN.out.perms      : null
            def scoring_accum_ch         = accum_results     ? accum_results.results               : null
            // genomic_info comes from params.gene_ensembl_file (resolved inside scoring.nf)
            // scoring_caas_* are built above, outside this block — see the comment there.

            // FOP per-pair PSS weights (contrast_hypotheses_pairs.tsv) for
            // domain-pooled scoring. Written by 4.Independent_contrasts.Rmd into
            // the Traitfiles dir and carried through CHECK_MIN_CONTRASTS's
            // traitfiles_ok_dir. When contrast selection ran, this channel is
            // authoritative: it emits the file, or NO_HYP_PAIRS for a
            // single-hypothesis run, and scoring.nf does not look in outdir.
            def scoring_hyp_pairs_ch = contrast_out
                ? (contrast_out.trait_dir_out ?: Channel.empty())
                      .map { d ->
                          def f = d ? file("${d}/contrast_hypotheses_pairs.tsv") : null
                          (f && f.exists()) ? f : file('NO_HYP_PAIRS')
                      }
                : null

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
                scoring_caas_pos_cycle_caas_ch, // per (gene,position,side,cycle) caas_score → p.emp
                scoring_caas_pos_sample_ch,  // cycle-stratified sample for report distribution plots
                scoring_caas_pos_quantiles_ch, // per (cycle,scheme) null distribution shape
                scoring_hyp_pairs_ch           // contrast_hypotheses_pairs.tsv — FOP domain-pool weights
            )
            ran_any = true

            // The evidence of the N best positions: what each domain of each hypothesis saw. It re-scores the rows of
            // those positions from the observed discovery.tab, so it needs the one this run wrote or was given.
            if (evidence_top_n > 0) {
                if (!evidence_inputs) {
                    error "caas_evidence_top_n > 0 needs the observed discovery.tab: run the alignments through the core (ct_tool 'discovery') or give --discovery_from."
                }
                CAAS_EVIDENCE(evidence_inputs.discovery.first(), SCORING.out.position_scores.first(),
                              evidence_inputs.design.first(), evidence_inputs.tree.first())
            }

            // CAAS_SIGNIFICANCE_REPORT: a DISTINCT, LATER stage than
            // CAAS_META_CAAS_REPORT (run above inside the run_meta_caas
            // block). It must run after SCORING because it joins
            // position_scores.tsv (p.emp/p.adj_bh/p.adj_sam) and gene_scores.tsv
            // (gene_caas_pperm/gene_caas_pperm_adj) onto the postproc-filtered
            // discovery table (the exact pooled dataset evaluated by SCORING),
            // with fallback to CT_META_CAAS's global_meta_caas.tsv when standalone.
            if (params.scoring) {
                def signif_caas_upstream = scoring_postproc_ch
                if (!signif_caas_upstream && core_observed) {
                    signif_caas_upstream = core_observed.global_meta_caas
                } else if (!signif_caas_upstream && observed_meta) {
                    signif_caas_upstream = observed_meta.global_meta_caas
                }

                if (signif_caas_upstream) {
                    CAAS_SIGNIFICANCE_REPORT(
                        signif_caas_upstream,
                        SCORING.out.position_scores,
                        SCORING.out.gene_scores
                    )
                    ran_any = true
                }
            }

            // VEP after SCORING: fed directly by position_scores.tsv.
            if (params.vep) {
                VEP(SCORING.out.position_scores)
            }

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
                def scoring_vep_pai_ch    = params.vep ? VEP.out.primateai_tsv : null
                def scoring_vep_cosmic_ch = params.vep ? VEP.out.cosmic_tsv    : null
                // POSENRICH background = caastools background.output (tested positions);
                // the engine restricts it to the cleaned_background genes.
                def posenrich_background_ch = core_observed ? core_observed.background : file('NO_FILE')

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
                    scoring_caas_pos_sample_ch,
                    scoring_caas_pos_cycle_caas_ch,
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

        // Standalone VEP is retired: VEP strictly depends on SCORING and is fed by position_scores.tsv.
        if (params.vep && !params.scoring) {
            error "VEP requires --scoring: standalone VEP has been retired. VEP depends on SCORING and is fed directly by position_scores.tsv."
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
