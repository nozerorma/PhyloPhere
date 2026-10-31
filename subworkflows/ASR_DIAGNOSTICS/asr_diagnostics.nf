// Report 9 describes completed observed scoring; it does not feed CT_POSTPROC.
def shellQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def diagnostics_source_hash = java.security.MessageDigest.getInstance('SHA-256').digest([
    'subworkflows/ASR_DIAGNOSTICS/local/9.ASR_diagnostics.Rmd',
    'subworkflows/ASR_DIAGNOSTICS/local/src/diagnostics_io.R',
    'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'
].collect { file("${baseDir}/${it}").text }.join('\n').getBytes('UTF-8')).encodeHex().toString()

process ASR_DIAGNOSTICS_REPORT {
    ext report_source_hash: diagnostics_source_hash
    cache 'deep'
    tag 'asr_diagnostics'
    label 'process_reporting'
    label 'error_retry'
    publishDir path: "${params.outdir}/asr_diagnostics", mode: 'copy', overwrite: true, pattern: '{tsv/**,plots/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    tuple val(observed_hash), path(disambiguation_dir, stageAs: 'inputs/observed/*')
    val posterior_threshold

    output:
    path '9.ASR_diagnostics.html', emit: report
    path 'tsv/**', emit: tables
    path 'plots/**', emit: plots, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ASR_DIAGNOSTICS/local"
    def style = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R"
    def report_params = groovy.json.JsonOutput.toJson([
        disambig_dir: disambiguation_dir.toString(),
        posterior_threshold: posterior_threshold, output_dir: '.'
    ])
    def rscript = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    echo '${task.ext.report_source_hash}' > report_source.sha256
    echo '${observed_hash}' > observed_input.sha256
    cp -R ${shellQuote(local_dir)}/. .
    cp ${shellQuote(style)} src/
    cat > report_params.json <<'PHYLOPHERE_REPORT_PARAMS'
${report_params}
PHYLOPHERE_REPORT_PARAMS
    REPORT_CORES=${task.cpus} ${rscript} - <<'PHYLOPHERE_RENDER_R'
    rmarkdown::render('9.ASR_diagnostics.Rmd',
        params = jsonlite::read_json('report_params.json', simplifyVector = TRUE),
        output_file = '9.ASR_diagnostics.html', envir = new.env())
PHYLOPHERE_RENDER_R
    """
}
