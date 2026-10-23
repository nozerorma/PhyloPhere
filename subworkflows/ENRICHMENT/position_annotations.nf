#!/usr/bin/env nextflow
// Annotation-only path: no positional test, permulation input or report 14.
process POSITION_ANNOTATIONS {
    label 'process_low'
    publishDir "${params.outdir}/posenrich/gmts", mode: 'copy', overwrite: true

    input:
    path gene_ensembl_file
    path domain_file
    path ucr_file
    path fubar_file
    path map_dir
    path background
    path fade_top
    path fade_bottom
    path custom_file

    output:
    path '*.gmt', emit: gmts
    path 'characterization_layers.tsv', emit: layers
    path 'position_characterization.tsv', emit: position_char

    script:
    def optional = { flag, source -> source.name.startsWith('NO_') ? '' : "--${flag} '${source}'" }
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/build_position_gmt.py \\
        ${optional('gene_ensembl_file', gene_ensembl_file)} \\
        ${optional('domain_variability_file', domain_file)} \\
        ${optional('ucr_positions_file', ucr_file)} \\
        ${optional('fubar_sites_file', fubar_file)} \\
        ${optional('map_dir', map_dir)} \\
        ${optional('cleaned_background', background)} \\
        ${optional('fade_sites_top_file', fade_top)} \\
        ${optional('fade_sites_bottom_file', fade_bottom)} \\
        ${optional('custom_marker_file', custom_file)} \\
        --output_dir .
    """
}
