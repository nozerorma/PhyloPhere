include { vepQuote } from './vep_contract'

process PRIMATEAI_MAP {
    cache 'deep'
    tag 'primateai'
    label 'process_long_compute'
    errorStrategy 'terminate'
    publishDir path: "${params.outdir}/vep", mode: 'copy', overwrite: true, pattern: 'primateai_mapped.tsv*'

    input:
    val fingerprint
    path targets
    path database
    path database_manifest

    output:
    path 'primateai_mapped.tsv', emit: primateai_tsv
    path 'primateai_mapped.tsv.manifest.json', emit: audit

    script:
    def src = "${baseDir}/subworkflows/VEP/local/src"
    """
    # Content fingerprint: ${fingerprint}
    cp ${vepQuote(src+'/map_to_primateai.py')} ${vepQuote(src+'/map_database.py')} ${vepQuote(src+'/vep_common.py')} .
    python3 map_to_primateai.py --targets ${vepQuote(targets)} --database ${vepQuote(database)} \
        --manifest ${vepQuote(database_manifest)} --output primateai_mapped.tsv
    """
}
