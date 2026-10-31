#!/usr/bin/env nextflow
// Report 8 consumes task artifacts directly; published directories are not inputs.
include { ppQuote; ppHash } from './ctpp_clustfilter'
def characterization_hash = ppHash([
    file("${baseDir}/subworkflows/CT_POSTPROC/local/8.Characterization_report.Rmd"),
    file("${baseDir}/subworkflows/CT_POSTPROC/local/src/report_io.R"),
    file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")])

process CT_POSTPROC_REPORT {
    ext source_hash: characterization_hash
    cache 'deep'
    tag 'caas_postproc_report'
    label 'process_reporting'
    publishDir path: "${params.outdir}/postproc", mode: 'copy', overwrite: true, pattern: '{CT_postproc_files/**,outliers/**,clusters/**,summary_statistics/**,disambiguation_characterization/**,postproc_inputs/**}'
    publishDir path: "${params.outdir}/html_reports", mode: 'copy', overwrite: true, pattern: '*.html'
    input:
    path discovery_file, stageAs: 'inputs/discovery/*'
    path filter_summary, stageAs: 'inputs/summary/*'
    path filter_files, stageAs: 'inputs/clusters/*'
    path selected_cluster, stageAs: 'inputs/selected/*'
    path gene_ensembl_file, stageAs: 'inputs/annotation/*'
    path gene_stats_file, stageAs: 'inputs/stats/*'
    path removed_genes, stageAs: 'inputs/removals/*'
    path filtered_discovery, stageAs: 'inputs/filtered/*'
    path quality_universe, stageAs: 'inputs/universe/*'
    output:
    path '8.Characterization_report.html', emit: report
    path 'outliers/**', emit: outliers
    path 'clusters/**', emit: clusters
    path 'summary_statistics/**', emit: summary_stats
    path 'disambiguation_characterization/**', emit: disambiguation_characterization
    path 'postproc_inputs/**', emit: gene_ensembl_input
    script:
    def local_dir = "${baseDir}/subworkflows/CT_POSTPROC/local"
    def settings = groovy.json.JsonOutput.toJson([
        processing_mode: params.caas_postproc_mode, gene_filter_mode: params.gene_filter_mode,
        gene_filter_strata: params.gene_filter_strata, gene_filter_k_far: params.gene_filter_k_far,
        gene_filter_k_near: params.gene_filter_k_near, remove_caas_clusters: params.remove_caas_clusters,
        selected_filter_minlen: params.filter_minlen, selected_filter_maxcaas: params.filter_maxcaas,
        selected_filter_min_positions: params.filter_min_positions,
        legacy_gene_count_gate: params.filter_min_positions == null,
        seed: params.seed ?: 1998, quality_universe: 'all_ct_candidates'])
    def report_params = groovy.json.JsonOutput.toJson([
        discovery_file: 'postproc_inputs/postproc_disambiguation_input.tsv',
        filtered_discovery_file: 'postproc_inputs/filtered_discovery.tsv',
        filter_summary_file: 'postproc_inputs/filter_summary.tsv', filter_dir: 'postproc_inputs/clusters',
        selected_cluster_file: 'postproc_inputs/selected_clusters.tsv',
        gene_ensembl_file: 'postproc_inputs/gene_annotation.tsv', gene_stats_file: 'postproc_inputs/gene_stats.tsv',
        removed_genes_file: 'postproc_inputs/removed_genes_summary.tsv',
        quality_universe_file: 'postproc_inputs/quality_universe.json',
        configuration_file: 'postproc_inputs/postproc_settings.json', output_dir: '.'])
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    def staged_clusters = filter_files instanceof Collection ? filter_files : [filter_files]
    def cluster_copies = staged_clusters.collect { "cp ${ppQuote(it)} postproc_inputs/clusters/" }.join('\n')
    """
    echo '${task.ext.source_hash}' > report_source.sha256
    cp -R ${ppQuote(local_dir)}/. .
    cp ${ppQuote("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")} src/
    mkdir -p postproc_inputs/clusters
    cp ${ppQuote(discovery_file)} postproc_inputs/postproc_disambiguation_input.tsv
    cp ${ppQuote(filtered_discovery)} postproc_inputs/filtered_discovery.tsv
    cp ${ppQuote(filter_summary)} postproc_inputs/filter_summary.tsv
    cp ${ppQuote(selected_cluster)} postproc_inputs/selected_clusters.tsv
    cp ${ppQuote(gene_ensembl_file)} postproc_inputs/gene_annotation.tsv
    cp ${ppQuote(gene_stats_file)} postproc_inputs/gene_stats.tsv
    cp ${ppQuote(removed_genes)} postproc_inputs/removed_genes_summary.tsv
    cp ${ppQuote(quality_universe)} postproc_inputs/quality_universe.json
    ${cluster_copies}
    cat > postproc_inputs/postproc_settings.json <<'PHYLOPHERE_SETTINGS'
${settings}
PHYLOPHERE_SETTINGS
    cat > report_params.json <<'PHYLOPHERE_PARAMS'
${report_params}
PHYLOPHERE_PARAMS
    REPORT_CORES=${task.cpus} ${run} - <<'PHYLOPHERE_RENDER'
    rmarkdown::render('8.Characterization_report.Rmd',
        params = jsonlite::read_json('report_params.json', simplifyVector = TRUE),
        output_file = '8.Characterization_report.html', envir = new.env())
PHYLOPHERE_RENDER
    """
}
