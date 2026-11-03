#!/usr/bin/env nextflow

// Descriptive position annotations, usable without the POSENRICH statistical workflow.
include { DOMAIN_VARIABILITY_GENERATION } from '../subworkflows/POSENRICH/domain_variability_generation.nf'
include { UCR_GENERATION } from '../subworkflows/POSENRICH/ucr_generation.nf'
include { POSITION_REGISTRY; POSENRICH_CHARACTERIZE } from '../subworkflows/POSENRICH/position_modules.nf'

def precomputedOrSentinel(String param_name, String param_value, String sentinel) {
    if (!param_value) return file(sentinel)
    if (file(param_value).exists()) return file(param_value)
    if (!param_value.endsWith('.gz') && file("${param_value}.gz").exists()) return file("${param_value}.gz")
    error "Supplied ${param_name} does not exist: ${param_value}"
}

// Descriptive annotations use SCORING inputs and do not need the enrichment universe or null blocks.
workflow POSITION_CHARACTERIZATION {
    take:
        position_scores
        postproc
        null_universe
        background_channel
        fade_sites_top_ch
        fade_sites_bottom_ch
    main:
        def scores_v = position_scores.collect().map { if (it.size()!=1) error('Expected exactly one SCORING position table'); it[0] }
        def observed_v = postproc.collect().map { if (it.size()!=1) error('Expected exactly one SCORING post-processing table'); it[0] }
        def background_v = (background_channel ?: Channel.empty()).ifEmpty { error('Position characterization requires a cleaned background (--scoring_background_input)') }.collect().map { if (it.size()!=1) error('Expected exactly one cleaned background'); it[0] }
        // Supplied tables take priority. Enabled layers with blank files are
        // generated once on the productive alignment columns.
        def domain_annotations_ch = Channel.value(file('NO_ANNOT_DOMAIN'))
        if (params.posenrich_domains) {
            if (params.domain_variability_file) {
                domain_annotations_ch = Channel.value(precomputedOrSentinel('domain_variability_file', params.domain_variability_file, 'NO_ANNOT_DOMAIN'))
            } else {
                DOMAIN_VARIABILITY_GENERATION(Channel.value(file(params.alignment)))
                domain_annotations_ch = DOMAIN_VARIABILITY_GENERATION.out.domain_variability_file
            }
        }
        def ucr_annotations_ch = Channel.value(file('NO_ANNOT_UCR'))
        if (params.posenrich_ucr) {
            if (params.ucr_positions_file) {
                ucr_annotations_ch = Channel.value(precomputedOrSentinel('ucr_positions_file', params.ucr_positions_file, 'NO_ANNOT_UCR'))
            } else {
                UCR_GENERATION(Channel.value(file(params.alignment)))
                ucr_annotations_ch = UCR_GENERATION.out.ucr_positions_file
            }
        }
        def fade_top_ch = (fade_sites_top_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_TOP') }
        def fade_bot_ch = (fade_sites_bottom_ch ?: Channel.empty()).ifEmpty { file('NO_FADE_BOTTOM') }

        def registry_out = POSITION_REGISTRY(
            scores_v,
            observed_v,
            (null_universe ?: Channel.empty()).ifEmpty { file('NO_NULL_UNIVERSE') },
            Channel.value(precomputedOrSentinel('caas_map_dir', params.caas_map_dir, 'NO_ANNOT_MAP')),
            background_v,
            Channel.value(precomputedOrSentinel('gene_ensembl_file', params.gene_ensembl_file, 'NO_ANNOT_ENSEMBL')),
            domain_annotations_ch,
            ucr_annotations_ch,
            Channel.value(precomputedOrSentinel('fubar_sites_file', params.fubar_sites_file, 'NO_ANNOT_FUBAR')),
            fade_top_ch, fade_bot_ch,
            Channel.value(precomputedOrSentinel('posenrich_custom_file', params.posenrich_custom_file, 'NO_ANNOT_CUSTOM')),
            Channel.value(file('NO_UNIVERSE')), Channel.value(file('NO_UNIVERSE_MANIFEST')),
            domain_annotations_ch.combine(ucr_annotations_ch)
                .combine(Channel.value(precomputedOrSentinel('fubar_sites_file', params.fubar_sites_file, 'NO_ANNOT_FUBAR')))
                .combine(Channel.value(precomputedOrSentinel('posenrich_custom_file', params.posenrich_custom_file, 'NO_ANNOT_CUSTOM')))
                .combine(fade_top_ch).combine(fade_bot_ch)
                .map { items -> items.withIndex().collect { p, i ->
                def m = file("${p}.manifest.json")
                m.exists() ? m : file("NO_SOURCE_MANIFEST_${i}")
                } }, Channel.value(file('NO_POSITION_CALIBRATION'))
        )
        def characterization_out = POSENRICH_CHARACTERIZE(registry_out.registry)

    emit:
        position_char = characterization_out.position_char
}

