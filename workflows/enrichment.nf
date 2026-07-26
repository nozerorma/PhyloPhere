#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: Functional Enrichment Workflow
 *
 * Runs gene-set and position-level functional enrichment analyses (FCS, STRING).
 */

include { SCORING_AMI_REPORT; SCORING_COMPARE_REPORT } from "${baseDir}/subworkflows/ENRICHMENT/scoring_enrichment.nf"
include { SCORING_FCS_REPORT }                            from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { RER_FCS_REPORT  as MODULE_FCS_RER }             from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { POSENRICH }                                      from "${baseDir}/subworkflows/ENRICHMENT/posenrich.nf"
include { DOMINO_MODULES }                                 from "${baseDir}/subworkflows/ENRICHMENT/domino.nf"

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
        background_output_ch
        vep_primateai_ch     // optional: PrimateAI-3D score TSV (null when --vep not run)
        vep_cosmic_ch        // optional: COSMIC Mutant Census score TSV (null when --vep not run)
        fade_sites_top_ch    // optional: fade_sites_top.csv (null when --fade not run)
        fade_sites_bottom_ch // optional: fade_sites_bottom.csv (null when --fade not run)

    main:
        def fcs_universe_ch = (cleaned_background_ch ?: Channel.empty())
            .ifEmpty { file('NO_BACKGROUND') }
            .collect()

        // RER has its own natural gene universe (every gene RERConverge actually
        // attempted), independent of the CAAS-discovery background. Falls back to
        // the CAAS background when no module-specific universe file is supplied,
        // preserving prior behaviour for runs that don't set this.
        def rer_universe_ch  = params.rer_universe_file  ? Channel.value(file(params.rer_universe_file))  : fcs_universe_ch

        def caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
            .ifEmpty { file(params.caas_perms_file ?: 'NO_FILE') }
            .collect()
            .map { it[0] }

        def caas_perm_scores_resolved = (caas_perm_scores_ch ?: Channel.empty())
            .ifEmpty { file('NO_FILE') }
            .collect()
            .map { it[0] }

        def fcs_out = SCORING_FCS_REPORT(fcs_stats, fcs_universe_ch, caas_perms_resolved)

        def annot_ch = fcs_stats.first()
        def trait_lbl = params.traitname ?: 'trait'

        def rer_perms_resolved = (rer_perms_ch ?: Channel.empty())
            .ifEmpty { file(params.rer_perms_file ?: 'NO_FILE') }
            .collect()
            .map { it[0] }

        def rer_fcs = MODULE_FCS_RER(
            Channel.value('scoring/rer'),   fcs_stats_rer,
            rer_universe_ch, Channel.value("13.FCS_rer_${trait_lbl}"),
            rer_perms_resolved, annot_ch)

        // AMI active module identification (DOMINO + STRING labelling)
        def run_ami = params.scoring_ami ?: params.scoring_string ?: false
        def ami_out = null
        if (run_ami) {
            def gene_lists_ch = gene_lists.ifEmpty { file('NO_GENE_LISTS') }
            def domino_out = DOMINO_MODULES(
                fcs_universe_ch,
                gene_lists_ch,
                params.domino_network_score_thr ?: 700,
                params.domino_slice_thr ?: 0.3,
                params.domino_module_thr ?: 0.05,
                params.string_db_dir ?: ''
            )
            ami_out = SCORING_AMI_REPORT(
                gene_lists_ch, fcs_universe_ch, gene_scores,
                domino_out.network_sif, domino_out.modules_dir
            )
        }

        // compare report
        def caas_all = fcs_out.fcs_all_results.ifEmpty { file('NO_FCS_ALL') }
        def rer_all  = rer_fcs.fcs_all_results.ifEmpty  { file('NO_FCS_ALL') }
        // Leading-edge tables feed the Comparison report's orthogonal composite
        // score (percentile concentration + cross-module corroboration).
        def caas_le  = fcs_out.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }
        def rer_le   = rer_fcs.fcs_leading_edge.ifEmpty { file('NO_LEADING_EDGE') }

        def cmp_out = SCORING_COMPARE_REPORT(caas_all, rer_all, caas_le, rer_le)

        def final_reports = fcs_out.report
            .mix(cmp_out.report)
            .mix(rer_fcs.report)
        if (run_string) {
            final_reports = final_reports.mix(string_out.report)
        }

        // ── Position-level FCS ──
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
            // Position background = caastools background.output (tested positions),
            // restricted downstream to the cleaned_background genes. An optional
            // param overrides the in-run channel for standalone re-runs.
            def pos_background_ch = params.posenrich_background_file ?
                Channel.fromPath(params.posenrich_background_file).ifEmpty { file('NO_FILE_POS_BG') } :
                (background_output_ch ? background_output_ch : file('NO_FILE_POS_BG'))

            // Position Characterisation (PrimateAI-3D + COSMIC validation, moved here
            // from the Scoring report) reuses gene_scores (already a take: param) plus
            // the same VEP outputs SCORING's own report renders with.
            def pos_gene_scores_ch = gene_scores ? gene_scores : file('NO_FILE_GENE_SCORES')
            def pos_vep_primateai_ch = (vep_primateai_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_VEP_PAI') }
            def pos_vep_cosmic_ch   = (vep_cosmic_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_VEP_COSMIC') }

            // FADE_top_sig / FADE_bottom_sig position group (the classic BF>=100
            // sites, always produced by FADE_JSON_TO_CSV when --fade runs).
            def pos_fade_sites_top_ch    = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_FADE_TOP') }
            def pos_fade_sites_bottom_ch = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FILE_FADE_BOTTOM') }

            POSENRICH(
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
            final_reports = final_reports.mix(POSENRICH.out.report)
        }

    emit:
        report = final_reports
}
