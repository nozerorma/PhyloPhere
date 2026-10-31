#!/usr/bin/env nextflow

/*
##
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: ct_postproc.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT Post-Processing Workflow: Handles cluster filtering and characterization of
 *  CT discovery results. Provides parameter sweep (exploratory mode) and single
 *  parameter filtering (filter mode) with optional report generation.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { CAAS_PREPARE_POSTPROC_INPUT; CT_FILTER; CT_FILTER_SUMMARY; CAAS_FILTER_GENES; CAAS_BACKGROUND_CLEANUP; clusterParameterGrid; clusterFileSuffix; ppHash } from '../subworkflows/CT_POSTPROC/ctpp_clustfilter'
include { CT_POSTPROC_REPORT } from '../subworkflows/CT_POSTPROC/ctpp_characterization'
include { ASR_DIAGNOSTICS } from './asr_diagnostics'

workflow CT_POSTPROC {
    take:
    disambiguation_input_channel
    background_genes_channel
    disambiguation_dir_channel
    hyp_pairs_channel
    taxid_map_channel   // curated tax_id map (NAME_CURATION) or null: then --tax_id is read as given

    main:
    assert disambiguation_input_channel != null || params.disambiguation_input : 'CT_POSTPROC requires observed master input'
    def discovery_ch = (disambiguation_input_channel ?: Channel.value(file(params.disambiguation_input, checkIfExists: true)))
        .collect().map { hits ->
            assert hits.size() == 1 : 'Expected exactly one observed master'
            hits[0]
        }
    def observed_ch = discovery_ch.map { master ->
        assert master.isFile() : 'CT_POSTPROC requires a master CSV'
        def stamp = master.parent.resolve('observed_scoring.json')
        tuple(ppHash([master, stamp]), master, stamp)
    }
    if (params.asr_diagnostics) {
        ASR_DIAGNOSTICS(disambiguation_dir_channel ?: discovery_ch.map { it.parent })
    }
    def pairs_ch
    if (hyp_pairs_channel != null) {
        pairs_ch = hyp_pairs_channel.map { p -> p.name.startsWith('NO_') ? [] : p }.ifEmpty([])
    } else {
        def p = file("${params.outdir}/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv")
        pairs_ch = Channel.value(p.exists() ? p : [])
    }
    def taxid_ch
    if (taxid_map_channel != null) {
        taxid_ch = taxid_map_channel.map { t -> t.name.startsWith('NO_') ? [] : t }.ifEmpty([])
    } else {
        taxid_ch = Channel.value(params.tax_id ? file(params.tax_id, checkIfExists: true) : [])
    }
    def alignment = params.alignment ? file(params.alignment, checkIfExists: true) : []
    def maps = params.caas_map_dir ? file(params.caas_map_dir, checkIfExists: true) : []
    def map_hash = ppHash(maps ? [maps] : [])
    def prepared = CAAS_PREPARE_POSTPROC_INPUT(observed_ch, pairs_ch, alignment,
        ppHash(alignment ? [alignment] : []), taxid_ch)
    def background_ch
    if (background_genes_channel != null) {
        background_ch = background_genes_channel
    } else {
        assert params.background_input : 'CT_POSTPROC requires --background_input or upstream background'
        def bg = file(params.background_input, checkIfExists: true)
        if (bg.isDirectory()) {
            def hits = files("${bg}/*background_genes*", type: 'file')
            assert hits.size() == 1 : "Expected exactly one background_genes file in ${bg}; found ${hits.size()}"
            bg = hits[0]
        }
        assert bg.isFile() : "Expected background file: ${bg}"
        background_ch = Channel.value(bg)
    }
    background_ch = background_ch.collect().map { hits ->
        assert hits.size() == 1 : 'Expected exactly one global background'
        hits[0]
    }
    def minlen = params.filter_minlen.toString().toInteger()
    def maxcaas = params.filter_maxcaas.toString().toDouble()
    assert minlen >= 1 && Double.isFinite(maxcaas) && maxcaas >= 0 && maxcaas <= 1 : 'Invalid cluster parameters'
    if (params.filter_min_positions != null) {
        assert params.filter_min_positions.toString() ==~ /[1-9][0-9]*/ : 'filter_min_positions must be a positive integer'
    }
    assert params.caas_postproc_mode in ['filter','exploratory'] : 'Invalid caas_postproc_mode'
    def combos = params.caas_postproc_mode == 'exploratory' ? clusterParameterGrid(
        params.minlen_values.split(',').collect { it.trim().toInteger() },
        params.maxcaas_values.split(',').collect { it.trim().toDouble() }, minlen, maxcaas) : [[minlen,maxcaas]]
    assert combos.every { it[0] >= 1 && Double.isFinite(it[1]) && it[1] >= 0 && it[1] <= 1 } : 'Invalid exploratory parameters'
    def combinations = Channel.fromList(combos).combine(prepared.prepared_discovery)
        .map { l,c,d -> tuple(params.caas_postproc_mode,l,c,d) }
    def filters = CT_FILTER(combinations, maps, map_hash)
    def all_filters = filters.filtered_files.collect()
    def summary = CT_FILTER_SUMMARY(all_filters)
    def selected = all_filters.map { files ->
        def hits = files.flatten().findAll { it.name.endsWith(clusterFileSuffix(minlen,maxcaas)) }
        assert hits.size() == 1 : 'Expected one selected cluster artifact'
        hits[0]
    }
    assert params.gene_ensembl_file : 'CT_POSTPROC requires --gene_ensembl_file'
    def lengths = file(params.gene_ensembl_file, checkIfExists: true)
    def selected_contract = filters.filter_contracts.collect().map { hits ->
        def found = hits.flatten().findAll { it.name.endsWith(clusterFileSuffix(minlen,maxcaas)+'.scoring.json') }
        assert found.size()==1 : 'Expected one selected cluster contract'
        found[0]
    }
    def genes = CAAS_FILTER_GENES(prepared.prepared_discovery, prepared.quality_universe, selected_contract, lengths, selected, maps, map_hash)
    def background = CAAS_BACKGROUND_CLEANUP(background_ch, genes.removed_genes)
    CT_POSTPROC_REPORT(prepared.prepared_discovery, summary.summary, all_filters, selected,
        lengths, genes.gene_stats, genes.removed_genes, genes.filtered_discovery, prepared.quality_universe)
    emit:
    filter_summary = summary.summary
    filter_dir = Channel.value("${params.outdir}/postproc")
    filtered_discovery = genes.filtered_discovery
    cleaned_background = background.cleaned_background_main
}
