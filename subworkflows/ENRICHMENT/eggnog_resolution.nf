#!/usr/bin/env nextflow
// eggnog_resolution.nf — Provide the eggNOG members and annotations files of POSENRICH when none is supplied.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  EGGNOG_RESOLUTION: provides the --egg_members_file / --egg_annotations_file pair when
 *  either is left blank (bin/resolve_eggnog.py). By default it copies the eggNOG 5.0
 *  Primates-level (taxid 9443) human-member pair versioned in subworkflows/ENRICHMENT/dat/.
 *  With params.auto_fetch_eggnog it downloads the pair of params.eggnog_taxid instead.
 *  Called from workflows/enrichment.nf.
 *
 *  Consumes:  no channel input (params.eggnog_taxid or params.clade_taxid, params.ref_species_taxid)
 *  Produces:  core_inputs/eggnog/ with the members and annotations files and eggnog_source.json,
 *             which records the origin and the checksums of the files used
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── eggNOG pair ──────────────────────────────────────────────────────────────

process RESOLVE_EGGNOG {
    tag "egg_members_file/egg_annotations_file"
    label 'process_low'

    publishDir "${params.outdir}/core_inputs/eggnog", mode: 'copy', overwrite: true

    output:
    path "*_members_*.tsv.gz", emit: egg_members_file
    path "*_annotations_*.tsv.gz", emit: egg_annotations_file
    path "eggnog_source.json", emit: source

    script:
    def tax_level = params.eggnog_taxid ?: (params.clade_taxid ?: "9443")
    def ref_taxid = params.ref_species_taxid ?: "9606"
    def fetch = params.auto_fetch_eggnog ? "--fetch" : ""
    """
    python3 ${baseDir}/bin/resolve_eggnog.py --output-dir . --tax-level "${tax_level}" --ref-taxid "${ref_taxid}" ${fetch}
    """

    stub:
    """
    echo '' | gzip > 9443_members_human.tsv.gz
    echo '' | gzip > 9443_annotations_human.tsv.gz
    echo '{}' > eggnog_source.json
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow EGGNOG_RESOLUTION {
    main:
        RESOLVE_EGGNOG()

    emit:
        egg_members_file = RESOLVE_EGGNOG.out.egg_members_file
        egg_annotations_file = RESOLVE_EGGNOG.out.egg_annotations_file
        source = RESOLVE_EGGNOG.out.source
}
