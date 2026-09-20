#!/usr/bin/env nextflow

/*
 * EGGNOG_RESOLUTION subworkflow
 *
 * Auto-generates POSENRICH's --egg_members_file / --egg_annotations_file pair
 * when either is left blank: fetches the eggNOG5 Primates-level (taxid 9443)
 * orthogroup members/annotations files, filtered to the subset
 * build_position_gmt.py actually reads (human "9606.ENSP*" members only —
 * every other member/annotation column is discarded on load anyway), falling
 * back to the vendored copies in assets/eggnog/ when offline
 * (bin/resolve_eggnog.py, same fetch-first/vendored-fallback pattern as
 * bin/resolve_gmts.py).
 */

process RESOLVE_EGGNOG {
    tag "auto-generate egg_members_file/egg_annotations_file"
    label 'process_low'

    publishDir "${params.outdir}/core_inputs/eggnog", mode: 'copy', overwrite: true

    output:
    path "9443_members_human.tsv.gz", emit: egg_members_file
    path "9443_annotations_human.tsv.gz", emit: egg_annotations_file

    stub:
    """
    echo '' | gzip > 9443_members_human.tsv.gz
    echo '' | gzip > 9443_annotations_human.tsv.gz
    """

    script:
    """
    python3 ${baseDir}/bin/resolve_eggnog.py --output-dir .
    """
}

workflow EGGNOG_RESOLUTION {
    main:
        RESOLVE_EGGNOG()

    emit:
        egg_members_file = RESOLVE_EGGNOG.out.egg_members_file
        egg_annotations_file = RESOLVE_EGGNOG.out.egg_annotations_file
}
