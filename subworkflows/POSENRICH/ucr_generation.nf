#!/usr/bin/env nextflow
// Global sliding UCR: the same Valdar/sliding implementation consumed by the
// positional collections. Explicit zero-based productive alignment coordinates.
process COMPUTE_UCR_POSITIONS {
    tag "auto-generate ucr_positions_file"
    label 'process_medium'
    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path alignment_dir

    output:
    path "ucr_positions.tsv.gz", emit: ucr_positions_file
    path "ucr_positions.tsv.gz.manifest.json", emit: manifest, optional: true

    script:
    """
    python3 ${baseDir}/bin/build_reference_tables.py alignment-ucr \\
        --alignment-dir "${alignment_dir}" \\
        --output ucr_positions.tsv.gz --cpus ${task.cpus}
    """

    stub:
    """
    printf 'gene\tposition\tregion_type\tmethod\tvariability\tcoordinate_space\n' | gzip > ucr_positions.tsv.gz
    """
}

workflow UCR_GENERATION {
    take:
        alignment_dir_ch
    main:
        COMPUTE_UCR_POSITIONS(alignment_dir_ch)
    emit:
        ucr_positions_file = COMPUTE_UCR_POSITIONS.out.ucr_positions_file
}
