#!/usr/bin/env nextflow
include { scoringHash } from '../subworkflows/SCORING/scoring_contract'

/*
 * PHYLOPHERE: CAAS Scoring Workflow
 *
 * Computes position-level and gene-level CAAS scores by integrating outputs
 * from CT_POSTPROC, FADE and RERConverge.
 * Runs once on the full postproc pool; directional characterisation is
 * performed post-scoring via the side column.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 * File: workflows/scoring.nf
 */

include { SCORING_EVIDENCE } from '../subworkflows/SCORING/scoring_evidence.nf'

include { SCORING_COMPUTE }                                                           from '../subworkflows/SCORING/scoring_compute.nf'
include { SCORING_NULL_POSITIONS }                                                    from '../subworkflows/SCORING/scoring_null_positions.nf'
include { SCORING_REPORT }                                                            from '../subworkflows/SCORING/scoring_report.nf'


workflow SCORING {

    take:
        postproc_ch              // Channel<path> or null — filtered_discovery.tsv
        fade_summary_top_ch      // Channel<path> or null — fade_summary_top.tsv
        fade_summary_bottom_ch   // Channel<path> or null — fade_summary_bottom.tsv
        rer_summary_ch           // Channel<path> or null — rerconverge_summary_{trait}.tsv
        genomic_info_ch          // Channel<path> or null — gene genomic coords TSV  (optional)
        fade_site_top_ch         // Channel<path> or null — fade_site_bf_top.tsv     (optional)
        fade_site_bot_ch         // Channel<path> or null — fade_site_bf_bottom.tsv  (optional)
        cleaned_background_ch    // Channel<path> or null — cleaned_background_main.txt (FCS universe)
        rer_perms_ch             // Channel<path> or null — RER permulation RDS (corStat) for RER FCS p.perm
        caas_perms_ch            // Channel<path> or null — optional CAAS compatibility export; FCS prefers blocks
        null_positions_ch        // Channel<path> or null, null_positions.parquet: the positions the null scores, with nd
        null_blocks_ch           // Channel<List<path>> or null, the blocks' directories (scored positions of the permuted cycles)
        null_counts_ch           // Channel<path> or null, counts.json of the null
        null_summary_ch          // Channel<tuple(cycle_summary.tsv, cycle_hist.tsv)> or null, the report's null panels
        caas_pos_sample_ch       // Channel<path> or null — cycle-stratified sample for report distribution plots
        caas_pos_quantiles_ch    // Channel<path> or null — per-cycle null distribution shape
        hypotheses_pairs_ch      // Channel<path> or null — contrast_hypotheses_pairs.tsv (FOP domain-pool weights)
        perm_manifest_ch         // Channel<path> or null — permulation_manifest.tsv (design of the canonical pairs of each null cycle)
        perm_harvest_ch          // Channel<path> or null — permulation_harvest.tsv (draws, rejections, acceptance by PSS tolerance, capacity)
        tested_positions_ch      // Channel<path> or null — caastools background.output (tested positions): limits the null to the positional universe
        cycle_roster_ch          // Channel<path> or null — cycle_roster.txt of the CAAS_CORE_MERGE run: N of the null, cycles with no detection included

    main:
        assert params.traitname : "SCORING requires --traitname"

        def optional_input = { value, sentinel -> value ? file(value,checkIfExists:true) : file(sentinel) }
        def score_dir = "${baseDir}/subworkflows/SCORING"
        def hash_sources = { names -> scoringHash(names.collect { file("${baseDir}/subworkflows/${it}") } + [file("${score_dir}/scoring_contract.nf")]) }
        def compute_hash = hash_sources(['SCORING/local/src/scoring_compute.R','SCORING/local/src/scoring_io.R','SCORING/local/src/scoring_input.py','SCORING/local/src/observed_core_scores.py','SCORING/local/src/aa_grouping.R','SCORING/local/src/null_position_p.R','CT_DISAMBIGUATION/local/src/core/scores.py'])
        def evidence_hash = hash_sources(['SCORING/local/src/scoring_evidence.R','SCORING/local/src/scoring_io.R','SCORING/local/src/caas_block_contract.R','ENRICHMENT/local/src/fcs_prepare_blocks.R','ENRICHMENT/local/src/gene_p.R','ENRICHMENT/local/src/caas_blocks.R'])
        def report_hash = hash_sources(['SCORING/local/10.Scoring_report.Rmd','SCORING/local/src/scoring_io.R','SCORING/local/src/report_evidence.R','TRAIT_ANALYSIS/local/src/report_style.R'])

        // ── Resolve each input with .ifEmpty{} fallbacks ───────────────────

        def resolved_postproc
        if (postproc_ch) {
            // Broadcastable value channel: filtered_discovery is consumed
            // by BOTH SCORING_COMPUTE and the report's functional-group overlay tab.
            resolved_postproc = postproc_ch.collect().map { assert it.size()==1 : "SCORING requires exactly one postproc table"; it[0] }
        } else if (params.scoring_postproc_input) {
            def f = file(params.scoring_postproc_input)
            assert f.exists() : "SCORING: postproc input not found: ${params.scoring_postproc_input}"
            resolved_postproc = Channel.value(f)
        } else {
            error("SCORING requires postproc input")
        }

        def resolved_fade_top = (fade_summary_top_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.scoring_fade_summary_top,'NO_FADE_TOP') }.collect().map { assert it.size()==1 : 'Expected one FADE top summary'; it[0] }

        def resolved_fade_bottom = (fade_summary_bottom_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.scoring_fade_summary_bottom,'NO_FADE_BOTTOM') }.collect().map { assert it.size()==1 : 'Expected one FADE bottom summary'; it[0] }

        def resolved_rer = (rer_summary_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.scoring_rer_input,'NO_RER') }.collect().map { assert it.size()==1 : 'Expected one RER summary'; it[0] }

        def resolved_genomic_info = (genomic_info_ch ?: Channel.empty())
            .ifEmpty {
                def gi = params.gene_ensembl_file ?: ''
                gi ? file(gi,checkIfExists:true) : file('NO_GENOMIC_INFO')
            }

        // .collect() turns these into value channels so they can be safely consumed
        // by both SCORING_COMPUTE and SCORING_REPORT below — without it, a regular
        // (queue) channel is drained by its first consumer, leaving the report with
        // nothing and silently dropping the FADE site data from the rendered report.
        // .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] } unwraps the single-element list .collect() produces back to
        // a scalar path, since both processes call `.name` directly on this input.
        def resolved_fade_site_top_ch = (fade_site_top_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.scoring_fade_site_top,'NO_FADE_SITE_TOP.txt') }
            .collect()
            .map { assert it.size()==1 : "Expected one FADE top site table"; it[0] }

        def resolved_fade_site_bot_ch = (fade_site_bot_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.scoring_fade_site_bottom,'NO_FADE_SITE_BOT.txt') }
            .collect()
            .map { assert it.size()==1 : "Expected one FADE bottom site table"; it[0] }

        def resolved_background = (cleaned_background_ch ?: Channel.empty())
            .collect()
            .ifEmpty {
                def bg = params.scoring_background_input ?: ''
                bg ? [file(bg,checkIfExists:true)] : [file('NO_BACKGROUND')]
            }
            .map { assert it.size()==1 : "Expected one cleaned background"; it[0] }

        // FOP per-pair PSS weights for domain-pooled scoring (contrast_hypotheses_pairs.tsv).
        // NO_HYP_PAIRS -> scoring_compute.R treats every position as single-hypothesis.
        // Resolution:
        //   1. --scoring_hypotheses_pairs;
        //   2. integrated run (hypotheses_pairs_ch given): the file emitted by this
        //      run's contrast selection, or NO_HYP_PAIRS for a single-hypothesis run;
        //   3. standalone scoring re-run: auto-discover in outdir, else NO_HYP_PAIRS.
        // An integrated run never reads outdir, so a file published by an earlier
        // run cannot stand in for this run's contrast selection.
        def hp_param = params.scoring_hypotheses_pairs ?: ''
        def resolved_hyp_pairs
        if (hp_param) {
            resolved_hyp_pairs = Channel.value(file(hp_param,checkIfExists:true))
        } else if (hypotheses_pairs_ch != null) {
            resolved_hyp_pairs = hypotheses_pairs_ch.ifEmpty(file('NO_HYP_PAIRS')).collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }
        } else {
            def auto = file("${params.outdir}/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv")
            resolved_hyp_pairs = Channel.value(auto.exists() ? auto : file('NO_HYP_PAIRS'))
        }

        // Design of the null cycles (permulation_manifest.tsv), for the design panel of the report.
        // NO_PERM_MANIFEST -> the report omits the panel. --caas_perm_manifest_file overrides.
        def pm_param = params.caas_perm_manifest_file ?: ''
        def resolved_perm_manifest
        if (pm_param) {
            resolved_perm_manifest = Channel.value(file(pm_param,checkIfExists:true))
        } else if (perm_manifest_ch != null) {
            resolved_perm_manifest = perm_manifest_ch.ifEmpty(file('NO_PERM_MANIFEST')).collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }
        } else {
            resolved_perm_manifest = Channel.value(file('NO_PERM_MANIFEST'))
        }

        // Audit of the null harvest (permulation_harvest.tsv), for the "Null harvest" tab of the report.
        // NO_PERM_HARVEST -> the report omits the tab. --caas_perm_harvest_file overrides.
        def ph_param = params.caas_perm_harvest_file ?: ''
        def resolved_perm_harvest
        if (ph_param) {
            resolved_perm_harvest = Channel.value(file(ph_param,checkIfExists:true))
        } else if (perm_harvest_ch != null) {
            resolved_perm_harvest = perm_harvest_ch.ifEmpty(file('NO_PERM_HARVEST')).collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }
        } else {
            resolved_perm_harvest = Channel.value(file('NO_PERM_HARVEST'))
        }

        // Optional compatibility export: resolve absence to a value so ENRICHMENT always runs.
        def caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
            .ifEmpty { file('NO_FILE') }
            .collect()
            .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }

        // Positional universe and N of the null (see null_position_stats.py). The tested positions are the live
        // background.output, else --posenrich_background_file (the same file POSENRICH Accumulation reads); the cleaned background
        // (resolved_background) limits the null to the genes of the run. Each falls back to a NO_* sentinel: the universe is
        // then every null-positive position, and the log says so.
        def resolved_tested_pos = (tested_positions_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.posenrich_background_file,'NO_TESTED_POSITIONS') }
            .collect()
            .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }

        // The position-level null: the live one of CAAS_CORE_MERGE_REDUCE, else the null of an earlier run
        // (--caas_null_dir: its caas_permulation/ directory). Without one every consumer receives the NO_* sentinels and
        // the position p-values stay NA.
        def null_dir = params.caas_null_dir ? file(params.caas_null_dir) : null
        if (null_dir) assert null_dir.isDirectory() && null_dir.resolve('null_positions.parquet').isFile() : 'SCORING requires a canonical --caas_null_dir'
        def prior_null = null_dir && file("${null_dir}/null_positions.parquet").exists()
        def have_null = (null_positions_ch != null) || prior_null
        def null_counters_v
        def null_fit_v
        def null_base_v
        def resolved_cycle_roster
        def null_summary_v
        if (have_null) {
            def null_positions_v = null_positions_ch != null ? null_positions_ch.collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] } : Channel.value(file("${null_dir}/null_positions.parquet"))
            // The block directories of a prior null. type: 'dir' is required: files() matches regular files by default
            // and would return an empty list here, so no blk_* is staged and the block validation fails downstream.
            def prior_blocks = null_blocks_ch == null ? files("${null_dir}/blocks/*", type: 'dir').sort() : []
            if (null_blocks_ch == null && !prior_blocks) error("--caas_null_dir has no block directories: ${null_dir}/blocks/")
            def null_blocks_v = null_blocks_ch != null ? null_blocks_ch : Channel.value(prior_blocks)
            def null_counts_v = null_counts_ch != null ? null_counts_ch.collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] } : Channel.value(file("${null_dir}/counts.json"))
            resolved_cycle_roster = cycle_roster_ch != null ? cycle_roster_ch.collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] } : Channel.value(file("${null_dir}/cycle_roster.txt"))
            // Preserve each summary/histogram pair; collect() otherwise flattens it
            // and the report's [0]/[1] indexing selects components of one Path.
            null_summary_v = null_summary_ch != null ? null_summary_ch.collect(flat: false).map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] } : Channel.value([file("${null_dir}/cycle_summary.tsv"), file("${null_dir}/cycle_hist.tsv")])
            def null_stats = SCORING_NULL_POSITIONS(null_positions_v, null_blocks_v, null_counts_v, resolved_postproc.collect().map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] },
                                                    resolved_tested_pos, resolved_background, resolved_cycle_roster)
            null_counters_v = null_stats.counters
            null_fit_v = null_stats.fit
            null_base_v = null_stats.base
        } else {
            null_counters_v = Channel.value(file('NO_NULL_COUNTERS'))
            null_fit_v = Channel.value(file('NO_NULL_FIT'))
            null_base_v = Channel.value(file('NO_NULL_BASE'))
            resolved_cycle_roster = Channel.value(file('NO_CYCLE_ROSTER'))
            null_summary_v = Channel.value([file('NO_NULL_SUMMARY'), file('NO_NULL_HIST')])
        }

        // ── Run scoring — single pass on full postproc pool ────────────────
        def compute_out = SCORING_COMPUTE(
            compute_hash,
            resolved_postproc,
            resolved_fade_top,
            resolved_fade_bottom,
            resolved_fade_site_top_ch,
            resolved_fade_site_bot_ch,
            resolved_rer,
            resolved_hyp_pairs,
            null_counters_v,
            null_fit_v,
            null_base_v,
            resolved_cycle_roster
        )

        // The same calibrated evidence feeds reporting, AMI and FCS, even with FCS disabled.
        // collect() emits nothing for an empty channel. Emit an empty list so the
        // saved-null fallback is evaluated and SCORING_EVIDENCE still runs without a null.
        def calibration_blocks = (null_blocks_ch ?: Channel.empty()).flatten().collect().ifEmpty { [] }.map { b ->
            b ? b : (params.caas_null_dir ? [file(params.caas_null_dir)] : [])
        }
        def evidence_out = SCORING_EVIDENCE(evidence_hash, compute_out.gene_scores, compute_out.position_scores,
            compute_out.fcs_stats, resolved_background, calibration_blocks, caas_perms_resolved)

        // ── Render report ──────────────────────────────────────────────────
        // Helper: for optional compute outputs, emit a sentinel when absent.
        def _opt = { ch, sentinel ->
            ch.ifEmpty { file(sentinel) }
        }

        def caas_pos_sample_resolved = (caas_pos_sample_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.caas_pos_sample_file,'NO_CAAS_POS_SAMPLE') }
            .collect()
            .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }
        def caas_pos_quantiles_resolved = (caas_pos_quantiles_ch ?: Channel.empty())
            .ifEmpty { optional_input(params.caas_pos_quantiles_file,'NO_CAAS_POS_QUANTILES') }
            .collect()
            .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }

        // N, universe size, unseen observed positions and saturation curve of the null (null_universe_summary.tsv);
        // the sentinel when no null was given, and the report then omits its section.
        def null_universe_resolved = compute_out.null_universe_summary
            .ifEmpty { file('NO_NULL_UNIVERSE') }
            .collect()
            .map { assert it.size()==1 : "SCORING expected exactly one artifact"; it[0] }

        // filtered_discovery (resolved_postproc, a value channel) = the observed
        // per-(gene,position) asr_path_score overlaid on the null;
        // when scoring runs it is always present (same condition as position_scores).
        def report_out = SCORING_REPORT(
            report_hash,
            evidence_out.position_scores,
            evidence_out.gene_scores,
            compute_out.gene_correlations,
            evidence_out.metadata,
            resolved_fade_site_top_ch,
            resolved_fade_site_bot_ch,
            resolved_genomic_info,
            null_summary_v.map { it[0] },
            null_summary_v.map { it[1] },
            caas_pos_sample_resolved,
            caas_pos_quantiles_resolved,
            resolved_postproc,
            resolved_background,
            resolved_perm_manifest,
            resolved_hyp_pairs,
            resolved_perm_harvest,
            null_universe_resolved
        )

        def final_reports = report_out.report

    emit:
        position_scores  = evidence_out.position_scores
        gene_scores      = evidence_out.gene_scores
        fcs_stats        = evidence_out.fcs_stats
        fcs_stats_rer    = compute_out.fcs_stats_rer
        fcs_stats_fade   = compute_out.fcs_stats_fade
        gene_lists       = evidence_out.gene_lists
        position_lists   = evidence_out.position_lists
        // WP10 POSENRICH modular: the registry consumes the evaluated positions and
        // the positional-null universe summary next to position_scores.tsv.
        postproc              = resolved_postproc
        null_universe_summary = compute_out.null_universe_summary
        report           = final_reports
        // The CAAS null used here; ENRICHMENT consumes this one, so the FCS p.perm and the scoring report share it.
        gene_evidence_metadata = evidence_out.metadata
        gene_evidence    = evidence_out.gene_evidence
        calibration      = evidence_out.calibration
        caas_perms       = caas_perms_resolved
}
