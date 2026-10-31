include { vepQuote } from './vep_contract'

process VEP_PREPARE {
    cache 'deep'
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir path: "${params.outdir}/vep/vep_inputs", mode: 'copy', overwrite: true

    input:
    val fingerprint
    path scores
    path maps
    path proteins

    output:
    path 'vep_targets.json', emit: targets
    path 'validated_inputs.json', emit: audit

    script:
    def py3 = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def src = "${baseDir}/subworkflows/VEP/local/src"
    def proteins_arg = proteins ? "--proteins ${vepQuote(proteins)}" : ''
    """
    cp ${vepQuote(src+'/vep_common.py')} ${vepQuote(src+'/prepare_targets.py')} .
    ${py3} prepare_targets.py --scores ${vepQuote(scores)} --maps ${vepQuote(maps)} \
        ${proteins_arg} --fingerprint ${vepQuote(fingerprint)}
    """
}
