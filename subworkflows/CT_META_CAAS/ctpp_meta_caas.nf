// CT_META_CAAS reports. Observed scoring owns meta_caas.tsv; these tasks
// publish the pattern HTML and canonical positive-score evidence HTML/TSVs.

def shellQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

// Expose these fingerprints through task.ext so Nextflow includes them in
// cache keys. Computing a digest only inside script: is insufficient.
def pattern_report_hash = java.security.MessageDigest.getInstance('SHA-256')
    .digest((file("${baseDir}/subworkflows/CT_META_CAAS/local/7.CAAS_pattern_annotation.Rmd").text + file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R").text).getBytes('UTF-8')).encodeHex().toString()
def significance_report_hash = java.security.MessageDigest.getInstance('SHA-256')
    .digest((file("${baseDir}/subworkflows/CT_META_CAAS/local/15.CAAS_significance_report.Rmd").text + file("${baseDir}/subworkflows/ENRICHMENT/local/src/report_presentation.R").text + file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R").text).getBytes('UTF-8')).encodeHex().toString()

process CAAS_META_CAAS_REPORT {
    ext report_source_hash: pattern_report_hash
    label 'process_reporting'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path discovery_input, stageAs: 'inputs/discovery/*'
    path background_input, stageAs: 'inputs/background/*'

    output:
    path '7.CAAS_pattern_annotation.html', emit: report

    script:
    def local_dir = "${baseDir}/subworkflows/CT_META_CAAS/local"
    def style = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R"
    def report_params = groovy.json.JsonOutput.toJson([
        discovery_input: discovery_input.toString(), background_input: background_input.toString(),
        traitname: (params.traitname ?: 'trait').toString()
    ])
    def rscript = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    echo '${task.ext.report_source_hash}' > report_source.sha256
    cp -R ${shellQuote(local_dir)}/. .
    mkdir -p src
    cp ${shellQuote(style)} src/
    cat > report_params.json <<'PHYLOPHERE_REPORT_PARAMS'
${report_params}
PHYLOPHERE_REPORT_PARAMS
    ${rscript} - <<'PHYLOPHERE_RENDER_R'
    rmarkdown::render('7.CAAS_pattern_annotation.Rmd',
        params = jsonlite::read_json('report_params.json', simplifyVector = TRUE),
        output_file = '7.CAAS_pattern_annotation.html', envir = new.env())
PHYLOPHERE_RENDER_R
    """
}

// Runs after SCORING. Producer p/q values and their original BH families are
// preserved. meta_caas remains an unused compatibility input for existing callers.
process CAAS_SIGNIFICANCE_REPORT {
    ext report_source_hash: significance_report_hash
    label 'process_reporting'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: 'significance/**'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path meta_caas, stageAs: 'inputs/meta/*'
    path position_scores, stageAs: 'inputs/positions/*'
    path gene_scores, stageAs: 'inputs/genes/*'
    path gene_metadata, stageAs: 'inputs/metadata/*'
    path gene_evidence, stageAs: 'inputs/evidence/*'
    path null_universe, stageAs: 'inputs/null/*'

    output:
    path '15.CAAS_significance_report.html', emit: report
    path 'significance/**', emit: significance_files
    path 'significance/meta_caas_significance.tsv', emit: meta_caas_significance

    script:
    def local_dir = "${baseDir}/subworkflows/CT_META_CAAS/local"
    def presentation = "${baseDir}/subworkflows/ENRICHMENT/local/src/report_presentation.R"
    def style = "${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R"
    def report_params = groovy.json.JsonOutput.toJson([
        meta_caas_input: meta_caas.toString(), position_scores_input: position_scores.toString(),
        gene_scores_input: gene_scores.toString(), gene_metadata_input: gene_metadata.toString(),
        gene_evidence_input: gene_evidence.toString(), null_universe_input: null_universe.toString(),
        traitname: (params.traitname ?: 'trait').toString(), seed: (params.seed != null ? params.seed : 1998).toString()
    ])
    def rscript = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    echo '${task.ext.report_source_hash}' > report_source.sha256
    cp -R ${shellQuote(local_dir)}/. .
    mkdir -p src
    cp ${shellQuote(presentation)} ${shellQuote(style)} src/
    cat > report_params.json <<'PHYLOPHERE_REPORT_PARAMS'
${report_params}
PHYLOPHERE_REPORT_PARAMS
    ${rscript} - <<'PHYLOPHERE_RENDER_R'
    rmarkdown::render('15.CAAS_significance_report.Rmd',
        params = jsonlite::read_json('report_params.json', simplifyVector = TRUE),
        output_file = '15.CAAS_significance_report.html', envir = new.env())
PHYLOPHERE_RENDER_R
    """
}
