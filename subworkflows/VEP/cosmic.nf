include { vepQuote } from './vep_contract'

process COSMIC_MAP {
    cache 'deep'
    tag 'cosmic'
    label 'process_long_compute'
    errorStrategy 'terminate'
    publishDir path: "${params.outdir}/vep", mode: 'copy', overwrite: true, pattern: 'cosmic_scores.tsv*'

    input:
    val fingerprint
    path targets
    path database
    path database_manifest

    output:
    path 'cosmic_scores.tsv', emit: cosmic_tsv
    path 'cosmic_scores.tsv.manifest.json', emit: audit

    script:
    def src = "${baseDir}/subworkflows/VEP/local/src"
    """
    # Content fingerprint: ${fingerprint}
    cp ${vepQuote(src+'/map_to_cosmic.py')} ${vepQuote(src+'/map_database.py')} ${vepQuote(src+'/vep_common.py')} .
    python3 map_to_cosmic.py --targets ${vepQuote(targets)} --database ${vepQuote(database)} \
        --manifest ${vepQuote(database_manifest)} --output cosmic_scores.tsv
    """
}
