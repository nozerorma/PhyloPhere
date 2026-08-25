#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: Functional Enrichment Workflow
 *
 * Runs gene-set and position-level functional enrichment analyses (FCS, AMI).
 */

include { SCORING_AMI_REPORT; SCORING_COMPARE_REPORT } from "${baseDir}/subworkflows/ENRICHMENT/scoring_enrichment.nf"
include { SCORING_FCS_REPORT }                            from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { RER_FCS_REPORT  as MODULE_FCS_RER }             from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { POSENRICH }                                      from "${baseDir}/subworkflows/ENRICHMENT/posenrich.nf"
// Nextflow forbids invoking the same process/subworkflow more than once in one workflow
// scope without a distinct alias per call site (DuplicateProcessInvocation) -- DOMINO_MODULES
// is called up to 3 times below (main/CAAS, FADE, RER), each against a different
// gene-list/universe/background, so each needs its own imported name.
include { DOMINO_MODULES as DOMINO_MODULES_MAIN }          from "${baseDir}/subworkflows/ENRICHMENT/domino.nf"
include { DOMINO_MODULES as DOMINO_MODULES_FADE }          from "${baseDir}/subworkflows/ENRICHMENT/domino.nf"
include { DOMINO_MODULES as DOMINO_MODULES_RER }           from "${baseDir}/subworkflows/ENRICHMENT/domino.nf"

// ── FADE Global background/gene-list union ──────────────────────────────────
// FADE's AMI network is built once per trait (not once per direction) against
// the union of both directions' own universes -- they differ by only ~12
// genes out of ~16k on real data, so the union is a safe superset either way.
// The Global gene list (for the unified 13.AMI_analysis.Rmd's "Global" FADE
// tab) is likewise the union of both directions' significant-gene lists.
process FADE_UNION_BACKGROUND {
    label 'process_low'

    input:
    path top_bg,    stageAs: 'fade_top_background.txt'
    path bottom_bg, stageAs: 'fade_bottom_background.txt'

    output:
    path 'fade_global_background.txt', emit: background

    script:
    """
    cat fade_top_background.txt fade_bottom_background.txt | sort -u > fade_global_background.txt
    """
}

process FADE_UNION_SIGNIFICANT {
    label 'process_low'

    input:
    path top_sig,    stageAs: 'fade_top_significant.txt'
    path bottom_sig, stageAs: 'fade_bottom_significant.txt'

    output:
    path 'fade_global_significant.txt', emit: significant

    script:
    """
    cat fade_top_significant.txt fade_bottom_significant.txt | sort -u > fade_global_significant.txt
    """
}

// ── Shared "universes" directory ─────────────────────────────────────────────
// Each tool (CAAS/FADE/RER) has its own correct gene universe -- they legitimately
// differ (verified on real data: RER ~16.1k, FADE ~16.1k-per-direction, CAAS's
// own cleaned_background a different count again) and none should be derived
// from another. Published together, trait-suffixed, purely for easy browsing/
// comparison -- plain copies, no new logic. NO_-prefixed sentinel inputs (tool
// didn't run) are simply skipped.
process PUBLISH_UNIVERSES {
    label 'process_low'
    publishDir path: "${params.outdir}/universes", mode: 'copy', overwrite: true

    input:
    path rer_bg
    path fade_bg
    path caas_bg
    val  traitname

    output:
    // Explicit dynamic filenames (not a "*.txt" glob) -- a glob here would also
    // match the Nextflow-staged INPUT files themselves (e.g. background.txt,
    // fade_global_background.txt, cleaned_background_main.txt all sit in the
    // task work dir as staged inputs), publishing confusing un-renamed
    // duplicates alongside the trait-suffixed copies below.
    path "rer_background_${traitname}.txt",             emit: rer_universe,  optional: true
    path "fade_background_${traitname}.txt",            emit: fade_universe, optional: true
    path "cleaned_background_main_${traitname}.txt",     emit: caas_universe, optional: true

    script:
    """
    rer_name=\$(basename ${rer_bg})
    fade_name=\$(basename ${fade_bg})
    caas_name=\$(basename ${caas_bg})
    [[ "\${rer_name}"  == NO_* ]]  || cp ${rer_bg}  "rer_background_${traitname}.txt"
    [[ "\${fade_name}" == NO_* ]]  || cp ${fade_bg} "fade_background_${traitname}.txt"
    [[ "\${caas_name}" == NO_* ]]  || cp ${caas_bg} "cleaned_background_main_${traitname}.txt"
    """
}

workflow ENRICHMENT {
    take:
        fcs_stats
        fcs_stats_rer
        fcs_stats_fade
        fcs_stats_accum
        gene_lists
        gene_scores
        cleaned_background_ch
        rer_perms_ch
        caas_perms_ch
        caas_perm_scores_ch
        caas_pos_pval_ch
        caas_pos_sample_ch
        position_scores
        position_lists       // SCORING's published position_lists/slice_{top,bottom,global}{25,10,5,1}.tsv dir
        background_output_ch
        vep_primateai_ch     // optional: PrimateAI-3D score TSV (null when --vep not run)
        vep_cosmic_ch        // optional: COSMIC Mutant Census score TSV (null when --vep not run)
        fade_sites_top_ch    // optional: fade_sites_top.csv (null when --fade not run)
        fade_sites_bottom_ch // optional: fade_sites_bottom.csv (null when --fade not run)
        rer_gene_lists_bg_ch        // RER's own background.txt (empty when RER didn't run)
        rer_gene_lists_interest_ch  // rer_significant/accelerating/decelerating.txt
        fade_gene_lists_bg_top_ch     // FADE top direction's own background.txt
        fade_gene_lists_bg_bottom_ch  // FADE bottom direction's own background.txt
        fade_gene_lists_sig_top_ch    // fade_top_significant.txt
        fade_gene_lists_sig_bottom_ch // fade_bottom_significant.txt

    main:
        def fcs_universe_ch = (cleaned_background_ch ?: Channel.empty())
            .ifEmpty { file('NO_BACKGROUND') }
            .collect()

        // rerconverge.nf/fade.nf emit these already `.collect()`-wrapped (each a
        // single-element List<Path>) since they're built via a `.filter(...).collect()`
        // chain there; normalize to a single Path here so every downstream use
        // (stageAs, DOMINO_MODULES' background arg, PUBLISH_UNIVERSES) sees one
        // consistent shape regardless of whether the tool ran this invocation.
        def unwrap1 = { ch -> ch.map { it instanceof List ? it[0] : it } }
        def rer_bg_r         = unwrap1(rer_gene_lists_bg_ch)
        def fade_bg_top_r    = unwrap1(fade_gene_lists_bg_top_ch)
        def fade_bg_bottom_r = unwrap1(fade_gene_lists_bg_bottom_ch)
        def fade_sig_top_r    = unwrap1(fade_gene_lists_sig_top_ch)
        def fade_sig_bottom_r = unwrap1(fade_gene_lists_sig_bottom_ch)

        // RER has its own natural gene universe (every gene RERConverge actually
        // attempted), independent of the CAAS-discovery background. Priority:
        // explicit --rer_universe_file override > RER's own live-run background
        // (rer_bg_r, from RER_MAIN.out.gene_lists_bg -- already computed above for
        // PUBLISH_UNIVERSES/DOMINO but previously never reached the FCS report) >
        // CAAS's background as a last-resort fallback for runs where RER didn't
        // run this invocation and no override was given (preserves prior behaviour).
        def rer_ran = (params.rer_tool || params.rer_continuous_file) as boolean
        def rer_universe_ch  = params.rer_universe_file
            ? Channel.value(file(params.rer_universe_file))
            : (rer_ran ? rer_bg_r : fcs_universe_ch)

        def caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
            .ifEmpty { file(params.caas_perms_file ?: 'NO_FILE') }
            .collect()
            .map { it[0] }

        def caas_perm_scores_resolved = (caas_perm_scores_ch ?: Channel.empty())
            .ifEmpty { file('NO_FILE') }
            .collect()
            .map { it[0] }

        // SCORING's own published gene_lists/ (scoring_compute.R) -- hoisted
        // here (rather than only inside the run_ami block below) since
        // SCORING_FCS_REPORT's gene_lists_dir param needs it too.
        def gene_lists_ch = gene_lists.ifEmpty { file('NO_GENE_LISTS') }

        def fcs_out = SCORING_FCS_REPORT(fcs_stats, fcs_universe_ch, caas_perms_resolved, gene_lists_ch)

        def annot_ch = fcs_stats
        def trait_lbl = params.traitname ?: 'trait'

        def rer_perms_resolved = (rer_perms_ch ?: Channel.empty())
            .ifEmpty { file(params.rer_perms_file ?: 'NO_FILE') }
            .collect()
            .map { it[0] }

        def rer_fcs = MODULE_FCS_RER(
            Channel.value('scoring/rer'),   fcs_stats_rer,
            rer_universe_ch, Channel.value("12.FCS_rer_${trait_lbl}"),
            rer_perms_resolved, annot_ch)

        // FADE and RER each have their own gene universe (verified to differ
        // meaningfully from CAAS's cleaned_background and from each other on
        // real data) -- computed once here regardless of whether AMI itself
        // runs, since PUBLISH_UNIVERSES below wants it too. NO_-prefixed
        // sentinels flow through when a tool didn't run this invocation.
        // FADE "ran" for AMI purposes whenever its gene lists/background are
        // actually available -- either a live --fade this invocation, or a
        // precomputed summary TSV supplied via --fade_json_dir_top/_bottom
        // (main.nf derives fade_gene_lists_bg/sig_* from that TSV the same way
        // it derives scoring_fade_top_ch/scoring_fade_bot_ch), mirroring how
        // rer_ran already treats rer_continuous_file as equivalent to rer_tool.
        def fade_ran = (params.fade || params.fade_json_dir_top || params.fade_json_dir_bottom) as boolean
        // rer_ran already hoisted above (needed earlier for rer_universe_ch).

        def fade_union_bg_ch
        def fade_union_sig_ch
        if (fade_ran) {
            def fade_union_bg  = FADE_UNION_BACKGROUND(fade_bg_top_r, fade_bg_bottom_r)
            def fade_union_sig = FADE_UNION_SIGNIFICANT(fade_sig_top_r, fade_sig_bottom_r)
            fade_union_bg_ch  = fade_union_bg.background
            fade_union_sig_ch = fade_union_sig.significant
        } else {
            fade_union_bg_ch  = Channel.value(file('NO_FADE_BACKGROUND'))
            fade_union_sig_ch = Channel.value(file('NO_FADE_SIGNIFICANT'))
        }
        // Manual override, mirroring rer_universe_ch above -- unlike RER (which
        // falls back to CAAS's background when unset), FADE's own union background
        // computed above is already correct on its own; this only matters when you
        // want to substitute a hand-picked universe file instead.
        if (params.fade_universe_file) {
            fade_union_bg_ch = Channel.value(file(params.fade_universe_file))
        }

        def rer_universe_for_publish_ch = rer_ran ? rer_bg_r : Channel.value(file('NO_RER_BACKGROUND'))
        PUBLISH_UNIVERSES(rer_universe_for_publish_ch, fade_union_bg_ch, fcs_universe_ch, Channel.value(trait_lbl))

        // AMI active module identification (DOMINO + STRING labelling)
        def run_ami = params.scoring_ami ?: params.scoring_string ?: false
        def ami_out = null
        if (run_ami) {
            def domino_out = DOMINO_MODULES_MAIN(
                fcs_universe_ch,
                gene_lists_ch,
                params.domino_network_score_thr ?: 700,
                params.domino_slice_thr ?: 0.3,
                params.domino_module_thr ?: 0.05,
                params.string_db_dir ?: '',
                'caas'
            )

            // FADE and RER each get their own DOMINO network, built against
            // the universe resolved above -- never reuse CAAS's network for
            // these. Each is gated on whether that tool actually ran this
            // invocation; when it didn't, a NO_-prefixed sentinel flows
            // through so 13.AMI_analysis.Rmd's NULL-arg checks (mirroring the
            // existing gene_scores/gs_arg pattern) skip that tool's section
            // gracefully.
            def fade_gene_list_files_arg
            def fade_background_arg
            def fade_domino_sif_arg
            def fade_domino_modules_arg
            def fade_domino_edges_arg
            if (fade_ran) {
                def fade_gene_lists_final_ch = fade_sig_top_r
                    .mix(fade_sig_bottom_r)
                    .mix(fade_union_sig_ch)
                    .collect()
                def fade_domino_out = DOMINO_MODULES_FADE(
                    fade_union_bg_ch,
                    fade_gene_lists_final_ch,
                    params.domino_network_score_thr ?: 700,
                    params.domino_slice_thr ?: 0.3,
                    params.domino_module_thr ?: 0.05,
                    params.string_db_dir ?: '',
                    'fade'
                )
                fade_gene_list_files_arg = fade_gene_lists_final_ch
                fade_background_arg      = fade_union_bg_ch
                fade_domino_sif_arg      = fade_domino_out.network_sif
                fade_domino_modules_arg  = fade_domino_out.modules_dir
                fade_domino_edges_arg    = fade_domino_out.edge_scores
            } else {
                fade_gene_list_files_arg = Channel.value(file('NO_FADE_GENE_LISTS'))
                fade_background_arg      = Channel.value(file('NO_FADE_BACKGROUND'))
                fade_domino_sif_arg      = Channel.value(file('NO_FADE_NETWORK_SIF'))
                fade_domino_modules_arg  = Channel.value(file('NO_FADE_MODULES_DIR'))
                fade_domino_edges_arg    = Channel.value(file('NO_FADE_EDGE_SCORES'))
            }

            def rer_gene_list_files_arg
            def rer_background_arg
            def rer_domino_sif_arg
            def rer_domino_modules_arg
            def rer_domino_edges_arg
            if (rer_ran) {
                def rer_domino_out = DOMINO_MODULES_RER(
                    rer_bg_r,
                    rer_gene_lists_interest_ch,
                    params.domino_network_score_thr ?: 700,
                    params.domino_slice_thr ?: 0.3,
                    params.domino_module_thr ?: 0.05,
                    params.string_db_dir ?: '',
                    'rer'
                )
                rer_gene_list_files_arg = rer_gene_lists_interest_ch
                rer_background_arg      = rer_bg_r
                rer_domino_sif_arg     = rer_domino_out.network_sif
                rer_domino_modules_arg = rer_domino_out.modules_dir
                rer_domino_edges_arg   = rer_domino_out.edge_scores
            } else {
                rer_gene_list_files_arg = Channel.value(file('NO_RER_GENE_LISTS'))
                rer_background_arg      = Channel.value(file('NO_RER_BACKGROUND'))
                rer_domino_sif_arg     = Channel.value(file('NO_RER_NETWORK_SIF'))
                rer_domino_modules_arg = Channel.value(file('NO_RER_MODULES_DIR'))
                rer_domino_edges_arg   = Channel.value(file('NO_RER_EDGE_SCORES'))
            }

            ami_out = SCORING_AMI_REPORT(
                gene_lists_ch, fcs_universe_ch, gene_scores,
                domino_out.network_sif, domino_out.modules_dir, domino_out.edge_scores,
                fade_gene_list_files_arg, fade_background_arg,
                fade_domino_sif_arg, fade_domino_modules_arg, fade_domino_edges_arg,
                rer_gene_list_files_arg, rer_background_arg,
                rer_domino_sif_arg, rer_domino_modules_arg, rer_domino_edges_arg,
                Channel.value(fade_ran), Channel.value(rer_ran)
            )
        }

        def final_reports = fcs_out.report
            .mix(rer_fcs.report)
        if (run_ami) {
            final_reports = final_reports.mix(ami_out.report)
        }

        // ── Position-level FCS ──
        def posenrich_out = null
        if (params.posenrich) {
            def pos_gene_ensembl_ch = Channel.fromPath(params.gene_ensembl_file).ifEmpty { file('NO_FILE') }
            def pos_domain_variability_ch = Channel.fromPath(params.domain_variability_file).ifEmpty { file('NO_FILE') }
            def pos_ucr_positions_ch = Channel.fromPath(params.ucr_positions_file).ifEmpty { file('NO_FILE') }
            def pos_fubar_sites_ch = Channel.fromPath(params.fubar_sites_file).ifEmpty { file('NO_FILE') }
            // Each of these is independently optional, and several can be absent
            // in the same run (e.g. no COSMIC/PrimateAI-3D DBs and no FADE run).
            // They must NOT share one literal 'NO_FILE' sentinel: Nextflow stages
            // path inputs by filename, so two absent inputs staged into the same
            // POSENRICH process call under the identical name 'NO_FILE' collide
            // ("input file name collision"). Every sentinel below is unique but
            // keeps the 'NO_FILE' prefix so the `=~ /^NO_FILE/` checks in
            // posenrich.nf still recognize it as absent.
            def pos_egg_members_ch = params.egg_members_file ? Channel.fromPath(params.egg_members_file).ifEmpty { file('NO_FILE_EGG_MEMBERS') } : file('NO_FILE_EGG_MEMBERS')
            def pos_egg_annotations_ch = params.egg_annotations_file ? Channel.fromPath(params.egg_annotations_file).ifEmpty { file('NO_FILE_EGG_ANNOT') } : file('NO_FILE_EGG_ANNOT')
            def pos_map_dir_ch = params.vep_map_dir ? Channel.fromPath(params.vep_map_dir).ifEmpty { file('NO_FILE_MAP_DIR') } : file('NO_FILE_MAP_DIR')
            def pos_cosmic_db_ch = params.cosmic_db ? Channel.fromPath(params.cosmic_db).ifEmpty { file('NO_FILE_COSMIC_DB') } : file('NO_FILE_COSMIC_DB')
            def pos_pai3d_db_ch = params.vep_primateai_db ? Channel.fromPath(params.vep_primateai_db).ifEmpty { file('NO_FILE_PAI3D_DB') } : file('NO_FILE_PAI3D_DB')
            def pos_cleaned_background_ch = cleaned_background_ch ? cleaned_background_ch : file('NO_FILE_BACKGROUND')

            def pos_caas_file_ch = position_scores ? position_scores : file('NO_FILE_CAAS')
            // SCORING's own published position percentile slices -- posenrich's
            // SOLE foreground source (no local re-ranking fallback: SCORING is a
            // mandatory upstream dependency of posenrich, not optional). If
            // position_lists is genuinely absent (e.g. --scoring off), this
            // sentinel flows through and posenrich_enrich.py hard-fails with a
            // clear message rather than silently degrading.
            def pos_lists_ch = (position_lists ?: Channel.empty()).ifEmpty { file('NO_FILE_POSITION_LISTS') }
            // Position background = caastools background.output (tested positions),
            // restricted downstream to the cleaned_background genes. An optional
            // param overrides the in-run channel for standalone re-runs.
            def pos_background_ch = params.posenrich_background_file ?
                Channel.fromPath(params.posenrich_background_file).ifEmpty { file('NO_FILE_POS_BG') } :
                (background_output_ch ? background_output_ch : file('NO_FILE_POS_BG'))

            // Position Characterisation (PrimateAI-3D + COSMIC validation, moved here
            // from the Scoring report) reuses gene_scores (already a take: param) plus
            // the same VEP outputs SCORING's own report renders with.
            // vep_primateai_ch/vep_cosmic_ch are null whenever --vep didn't run
            // this invocation (main.nf: `params.vep ? VEP.out.primateai_tsv :
            // null`) -- fall back to the same --scoring_vep_primateai/
            // --scoring_vep_cosmic precomputed-file params SCORING's own
            // report already resolves (workflows/scoring.nf's
            // resolved_vep_primateai/resolved_vep_cosmic), so a precomputed-
            // VEP run (--scoring_vep_primateai set, no live --vep) still
            // reaches Position Characterisation instead of silently going NULL.
            def pos_gene_scores_ch = gene_scores ? gene_scores : file('NO_FILE_GENE_SCORES')
            def pos_vep_primateai_ch = (vep_primateai_ch ?: Channel.empty())
                .ifEmpty { params.scoring_vep_primateai ? file(params.scoring_vep_primateai) : file('NO_FILE_VEP_PAI') }
            def pos_vep_cosmic_ch   = (vep_cosmic_ch ?: Channel.empty())
                .ifEmpty { params.scoring_vep_cosmic ? file(params.scoring_vep_cosmic) : file('NO_FILE_VEP_COSMIC') }

            // FADE_top_sig / FADE_bottom_sig position group (the classic BF>=100
            // sites, always produced by FADE_JSON_TO_CSV when --fade runs).
            def pos_fade_sites_top_ch    = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_FADE_TOP') }
            def pos_fade_sites_bottom_ch = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_FADE_BOTTOM') }

            posenrich_out = POSENRICH(
                pos_gene_ensembl_ch,
                pos_domain_variability_ch,
                pos_ucr_positions_ch,
                pos_fubar_sites_ch,
                pos_egg_members_ch,
                pos_egg_annotations_ch,
                pos_map_dir_ch,
                pos_cosmic_db_ch,
                pos_pai3d_db_ch,
                pos_cleaned_background_ch,
                pos_caas_file_ch,
                pos_lists_ch,
                pos_background_ch,
                annot_ch,
                params.posenrich_min_size,
                params.posenrich_max_size,
                pos_gene_scores_ch,
                pos_vep_primateai_ch,
                pos_vep_cosmic_ch,
                pos_gene_ensembl_ch,
                pos_fade_sites_top_ch,
                pos_fade_sites_bottom_ch
            )
            final_reports = final_reports.mix(posenrich_out.report)
        }

        // ── Comparison report — pulls in CAAS/RER FCS, plus (when available) AMI's
        // module descriptions and posenrich's position-level results, all called
        // above so their outputs exist as real objects by this point in the DAG. ──
        def caas_all = fcs_out.fcs_all_results.ifEmpty { file('NO_FCS_ALL') }
        def rer_all  = rer_fcs.fcs_all_results.ifEmpty  { file('NO_FCS_ALL') }
        // Leading-edge tables feed the Comparison report's Per-module
        // significance / Interesting Genes-Positions tables.
        def caas_le  = fcs_out.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }
        def rer_le   = rer_fcs.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }
        def caas_le_comp = fcs_out.fcs_leading_edge_composition.ifEmpty { file('NO_LEADING_EDGE_COMPOSITION') }
        def rer_le_comp  = rer_fcs.fcs_leading_edge_composition.ifEmpty  { file('NO_LEADING_EDGE_COMPOSITION') }

        def ami_module_desc_ch = (run_ami ? ami_out.module_descriptions : Channel.empty())
            .ifEmpty { file('NO_AMI_MODULE_DESC') }
        def ami_term_membership_ch = (run_ami ? ami_out.term_membership : Channel.empty())
            .ifEmpty { file('NO_AMI_TERM_MEMBERSHIP') }
        def posenrich_dotplot_ch = (posenrich_out ? posenrich_out.overall_dotplot : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_DOTPLOT') }
        def posenrich_le_ch = (posenrich_out ? posenrich_out.leading_edge : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_LE') }
        def posenrich_le_summary_ch = (posenrich_out ? posenrich_out.leading_edge_summary : Channel.empty())
            .ifEmpty { file('NO_POSENRICH_LE_SUMMARY') }
        // Interesting Genes/Positions tables reuse the SAME optional gene/position-
        // level inputs POSENRICH's own report already consumes (fcs_stats,
        // position_scores, VEP, FADE sites) -- all workflow-level take: params, so
        // available here whether or not POSENRICH itself ran this invocation.
        // Reuses annot_ch (fcs_stats.first(), already established above as the
        // broadcastable form of fcs_stats -- see its own definition) rather than
        // touching the raw fcs_stats channel a further time.
        def cmp_fcs_stats_ch  = annot_ch ?: file('NO_FCS_STATS')
        def cmp_pos_scores_ch = position_scores ? position_scores : file('NO_POS_SCORES')
        // Same precomputed-VEP fallback as pos_vep_primateai_ch/pos_vep_cosmic_ch above.
        def cmp_vep_pai_ch    = (vep_primateai_ch ?: Channel.empty())
            .ifEmpty { params.scoring_vep_primateai ? file(params.scoring_vep_primateai) : file('NO_VEP_PAI') }
        def cmp_vep_cosmic_ch = (vep_cosmic_ch ?: Channel.empty())
            .ifEmpty { params.scoring_vep_cosmic ? file(params.scoring_vep_cosmic) : file('NO_VEP_COSMIC') }
        def cmp_fade_top_ch   = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_TOP') }
        def cmp_fade_bot_ch   = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_BOTTOM') }
        // SCORING's published percentile slices -- gene_lists_ch already
        // established above (hoisted for SCORING_FCS_REPORT); position_lists
        // has no other consumer in this workflow yet, so it's resolved here.
        def cmp_position_lists_ch = position_lists ? position_lists.ifEmpty { file('NO_POSITION_LISTS') } : file('NO_POSITION_LISTS')

        def cmp_out = SCORING_COMPARE_REPORT(caas_all, rer_all, caas_le, rer_le,
                                              caas_le_comp, rer_le_comp,
                                              ami_module_desc_ch, ami_term_membership_ch,
                                              posenrich_dotplot_ch, posenrich_le_ch,
                                              posenrich_le_summary_ch,
                                              cmp_fcs_stats_ch, cmp_pos_scores_ch,
                                              cmp_vep_pai_ch, cmp_vep_cosmic_ch,
                                              cmp_fade_top_ch, cmp_fade_bot_ch,
                                              gene_lists_ch, cmp_position_lists_ch)
        final_reports = final_reports.mix(cmp_out.report)

    emit:
        report = final_reports
}
