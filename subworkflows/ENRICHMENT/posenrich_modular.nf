#!/usr/bin/env nextflow
// POSENRICH modular (WP10): common position registry + characterization +
// within-protein localization + global association + proximity.
// Gated by params.posenrich_modular; the historical POSENRICH inference stays omitted.

process POSITION_REGISTRY {
    label 'process_low'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'registry'

    input:
    path position_scores
    path postproc
    path null_universe
    path map_dir
    path background
    path gene_ensembl
    path domain_file
    path ucr_file
    path fubar_file
    path fade_top
    path fade_bottom
    path custom_file

    output:
    path 'registry', emit: registry

    script:
    // Same optional-input contract as POSITION_ANNOTATIONS: NO_* sentinels and
    // dangling paths mean absent; the registry tolerates every absent source.
    def optional = { flag, source -> source.name.startsWith('NO_') ? '' : "--${flag} '${source}'" }
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_registry.py \\
        ${optional('position_scores', position_scores)} \\
        ${optional('postproc', postproc)} \\
        ${optional('null_universe', null_universe)} \\
        ${optional('map_dir', map_dir)} \\
        ${optional('cleaned_background', background)} \\
        ${optional('gene_ensembl_file', gene_ensembl)} \\
        ${optional('domain_variability_file', domain_file)} \\
        ${optional('ucr_positions_file', ucr_file)} \\
        ${optional('fubar_sites_file', fubar_file)} \\
        ${optional('fade_sites_top_file', fade_top)} \\
        ${optional('fade_sites_bottom_file', fade_bottom)} \\
        ${optional('custom_marker_file', custom_file)} \\
        --output_dir registry --gmt_dir registry/gmts
    """
}

process POSENRICH_CHARACTERIZE {
    label 'process_low'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'characterization'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'characterization', emit: characterization

    script:
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_characterize.py \\
        --registry_dir registry --output_dir characterization
    """
}

process POSENRICH_LOCALIZATION {
    label 'process_low'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'localization'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'localization', emit: localization

    script:
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_localization.py \\
        --registry_dir registry --output_dir localization \\
        --ranking '${params.posenrich_ranking ?: 'caas_score'}' \\
        --aggregation '${params.posenrich_aggregation ?: 'mean'}' \\
        --n_permutations ${params.posenrich_n_permutations ?: 999} \\
        --seed ${params.posenrich_seed ?: 1}
    """
}

process POSENRICH_GLOBAL {
    label 'process_low'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'global'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'global', emit: global

    script:
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_global.py \\
        --registry_dir registry --output_dir global \\
        --ranking '${params.posenrich_ranking ?: 'caas_score'}'
    """
}

process POSENRICH_PROXIMITY {
    label 'process_low'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: 'proximity'

    input:
    path registry, stageAs: 'registry'

    output:
    path 'proximity', emit: proximity

    script:
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/posenrich_proximity.py \\
        --registry_dir registry --output_dir proximity \\
        --radios '${params.posenrich_radios ?: '0,5,10'}'
    """
}

process POSENRICH_MODULAR_REPORT {
    tag "posenrich_modular|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'
    publishDir path: "${params.outdir}/posenrich", mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path registry,       stageAs: 'registry'
    path characterization, stageAs: 'characterization'
    path localization,   stageAs: 'localization'
    path global_assoc,   stageAs: 'global'
    path proximity,      stageAs: 'proximity'

    output:
    path "POSENRICH_modular_report_${params.traitname ?: 'unknown_trait'}.html", emit: report

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    """
    cp -R ${local_dir}/* .
    Rscript -e "
        rmarkdown::render(
            'posenrich_modular_report.Rmd',
            params = list(
                registry_dir         = 'registry',
                characterization_dir = 'characterization',
                localization_dir     = 'localization',
                global_dir           = 'global',
                proximity_dir        = 'proximity',
                traitname            = '${traitname}',
                ranking              = '${params.posenrich_ranking ?: 'caas_score'}',
                aggregation          = '${params.posenrich_aggregation ?: 'mean'}',
                n_permutations       = ${params.posenrich_n_permutations ?: 999},
                seed                 = ${params.posenrich_seed ?: 1},
                radios               = '${params.posenrich_radios ?: '0,5,10'}'
            ),
            output_file = 'POSENRICH_modular_report_${traitname}.html'
        )
    "
    """
}

workflow POSENRICH_MODULAR {
    take:
    position_scores
    postproc
    null_universe
    map_dir
    background
    gene_ensembl
    domain_file
    ucr_file
    fubar_file
    fade_top
    fade_bottom
    custom_file

    main:
    def registry_out = POSITION_REGISTRY(position_scores, postproc, null_universe, map_dir,
                                        background, gene_ensembl, domain_file, ucr_file,
                                        fubar_file, fade_top, fade_bottom, custom_file)
    def characterization_out = POSENRICH_CHARACTERIZE(registry_out.registry)
    def localization_out   = POSENRICH_LOCALIZATION(registry_out.registry)
    def global_out         = POSENRICH_GLOBAL(registry_out.registry)
    def proximity_out      = POSENRICH_PROXIMITY(registry_out.registry)
    def report_out = POSENRICH_MODULAR_REPORT(registry_out.registry, characterization_out.characterization,
                                              localization_out.localization, global_out.global,
                                              proximity_out.proximity)

    emit:
    registry         = registry_out.registry
    characterization = characterization_out.characterization
    localization     = localization_out.localization
    global           = global_out.global
    proximity        = proximity_out.proximity
    report           = report_out.report
}
