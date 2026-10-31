include { vepQuote } from './vep_contract'

process ENSEMBL_VEP_ANNOTATE {
    // A live MySQL database is not a content-addressable snapshot.
    cache 'deep'
    tag 'ensembl_vep'
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir path: "${params.outdir}/vep", mode: 'copy', overwrite: true, pattern: 'ensembl_vep_mapped.tsv*'
    publishDir path: "${params.outdir}/vep/ensembl_raw", mode: 'copy', overwrite: true,
        pattern: '{candidate_variants.vcf,vep_tab_output.txt,vep_version.txt,vep_*.log,vep_warnings.txt}'

    input:
    val fingerprint
    path targets
    path cache_dir
    val settings

    output:
    path 'ensembl_vep_mapped.tsv', emit: ensembl_vep_tsv
    path 'ensembl_vep_mapped.tsv.manifest.json', emit: audit
    path 'candidate_variants.vcf', emit: variants
    path 'vep_tab_output.txt', optional: true, emit: raw
    path 'vep_version.txt', optional: true
    path 'vep_*.log', optional: true
    path 'vep_warnings.txt', optional: true

    script:
    def src = "${baseDir}/subworkflows/VEP/local/src"
    def cache_arg = cache_dir ? "--cache ${vepQuote(cache_dir)}" : ''
    """
    cp ${vepQuote(src+'/run_ensembl.py')} ${vepQuote(src+'/build_vep_variants.py')} \
        ${vepQuote(src+'/join_vep_output.py')} ${vepQuote(src+'/vep_common.py')} .
    python3 run_ensembl.py --targets ${vepQuote(targets)} --output ensembl_vep_mapped.tsv \
        --fingerprint ${vepQuote(fingerprint)} --mode ${vepQuote(settings.mode)} ${cache_arg} \
        --host ${vepQuote(settings.host)} --port ${vepQuote(settings.port)}
    """
}
