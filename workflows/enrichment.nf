#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: Functional Enrichment Workflow
 *
 * Runs gene-set and position-level functional enrichment analyses (FCS, ORA, STRING).
 */

include { SCORING_STRING_REPORT; SCORING_COMPARE_REPORT } from "${baseDir}/subworkflows/ENRICHMENT/scoring_enrichment.nf"
include { SCORING_FCS_REPORT }                            from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { RER_FCS_REPORT  as MODULE_FCS_RER }             from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { TOOL_FCS_REPORT as MODULE_FCS_FADE;
          TOOL_FCS_REPORT as MODULE_FCS_ACCUM }           from "${baseDir}/subworkflows/ENRICHMENT/fcs.nf"
include { POSENRICH }                                      from "${baseDir}/subworkflows/ENRICHMENT/posenrich.nf"

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

    main:
        def fcs_universe_ch = (cleaned_background_ch ?: Channel.empty())
            .ifEmpty { file('NO_BACKGROUND') }
            .collect()

        def caas_perms_resolved = (caas_perms_ch ?: Channel.empty())
            .collect()
            .map { it && it.size() > 0 ? it[0] : file(params.caas_perms_file ?: 'NO_FILE') }

        def caas_perm_scores_resolved = (caas_perm_scores_ch ?: Channel.empty())
            .collect()
            .map { it && it.size() > 0 ? it[0] : file('NO_FILE') }

        def fcs_out = SCORING_FCS_REPORT(fcs_stats, fcs_universe_ch, caas_perms_resolved)

        def annot_ch = fcs_stats.first()
        def trait_lbl = params.traitname ?: 'trait'

        def rer_perms_resolved = (rer_perms_ch ?: Channel.empty())
            .collect()
            .map { it && it.size() > 0 ? it[0] : file(params.rer_perms_file ?: 'NO_FILE') }

        def rer_fcs = MODULE_FCS_RER(
            Channel.value('scoring/rer'),   fcs_stats_rer,
            fcs_universe_ch, Channel.value("13.FCS_rer_${trait_lbl}"),
            rer_perms_resolved, annot_ch)

        def fade_fcs = MODULE_FCS_FADE(
            Channel.value('scoring/fade'),  fcs_stats_fade,
            fcs_universe_ch, Channel.value("13.FCS_fade_${trait_lbl}"),  annot_ch)

        def accum_fcs = MODULE_FCS_ACCUM(
            Channel.value('scoring/accum'), fcs_stats_accum,
            fcs_universe_ch, Channel.value("13.FCS_accum_${trait_lbl}"), annot_ch)

        // STRING ppi network
        def run_string = params.scoring_string ?: false
        def string_out = null
        if (run_string) {
            def gene_lists_ch = gene_lists.ifEmpty { file('NO_GENE_LISTS') }
            string_out = SCORING_STRING_REPORT(gene_lists_ch, fcs_universe_ch, gene_scores)
        }

        // compare report
        def caas_all  = fcs_out.fcs_all_results.ifEmpty   { file('NO_FCS_ALL') }
        def rer_all   = rer_fcs.fcs_all_results.ifEmpty   { file('NO_FCS_ALL') }
        def fade_all  = fade_fcs.fcs_all_results.ifEmpty  { file('NO_FCS_ALL') }
        def accum_all = accum_fcs.fcs_all_results.ifEmpty { file('NO_FCS_ALL') }

        def string_tsv_ch = run_string
            ? string_out.string_enrichment_tsvs
                .ifEmpty { file('NO_STRING_TSV') }
                .collect()
            : Channel.value([ file('NO_STRING_TSV') ])

        def cmp_out = SCORING_COMPARE_REPORT(caas_all, rer_all, fade_all, accum_all, string_tsv_ch)

        def final_reports = fcs_out.report
            .mix(cmp_out.report)
            .mix(rer_fcs.report)
            .mix(fade_fcs.report)
            .mix(accum_fcs.report)
        if (run_string) {
            final_reports = final_reports.mix(string_out.report)
        }

        // ── Position-level FCS ──
        if (params.posenrich) {
            def pos_gene_ensembl_ch = Channel.fromPath(params.gene_ensembl_file).ifEmpty { file('NO_FILE') }
            def pos_domain_variability_ch = Channel.fromPath(params.domain_variability_file).ifEmpty { file('NO_FILE') }
            def pos_ucr_positions_ch = Channel.fromPath(params.ucr_positions_file).ifEmpty { file('NO_FILE') }
            def pos_fubar_sites_ch = Channel.fromPath(params.fubar_sites_file).ifEmpty { file('NO_FILE') }
            def pos_egg_members_ch = params.egg_members_file ? Channel.fromPath(params.egg_members_file).ifEmpty { file('NO_FILE') } : file('NO_FILE')
            def pos_egg_annotations_ch = params.egg_annotations_file ? Channel.fromPath(params.egg_annotations_file).ifEmpty { file('NO_FILE') } : file('NO_FILE')
            def pos_map_dir_ch = params.vep_map_dir ? Channel.fromPath(params.vep_map_dir).ifEmpty { file('NO_FILE') } : file('NO_FILE')
            def pos_cosmic_db_ch = params.cosmic_db ? Channel.fromPath(params.cosmic_db).ifEmpty { file('NO_FILE') } : file('NO_FILE')
            def pos_cleaned_background_ch = cleaned_background_ch ? cleaned_background_ch : file('NO_FILE')

            def pos_caas_file_ch = position_scores ? position_scores : file('NO_FILE')
            // Position background = caastools background.output (tested positions),
            // restricted downstream to the cleaned_background genes. An optional
            // param overrides the in-run channel for standalone re-runs.
            def pos_background_ch = params.posenrich_background_file ?
                Channel.fromPath(params.posenrich_background_file).ifEmpty { file('NO_FILE') } :
                (background_output_ch ? background_output_ch : file('NO_FILE'))

            POSENRICH(
                pos_gene_ensembl_ch,
                pos_domain_variability_ch,
                pos_ucr_positions_ch,
                pos_fubar_sites_ch,
                pos_egg_members_ch,
                pos_egg_annotations_ch,
                pos_map_dir_ch,
                pos_cosmic_db_ch,
                pos_cleaned_background_ch,
                pos_caas_file_ch,
                pos_background_ch,
                annot_ch,
                params.posenrich_min_size,
                params.posenrich_max_size
            )
            final_reports = final_reports.mix(POSENRICH.out.report)
        }

    emit:
        report = final_reports
}
