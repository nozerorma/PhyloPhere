#!/usr/bin/env nextflow
include { RER_TRAIT } from '../subworkflows/RERCONVERGE/rer_trait'
include { RER_TREES } from '../subworkflows/RERCONVERGE/rer_trees'
include { RER_MATRIX } from '../subworkflows/RERCONVERGE/rer_matrix'
include { RER_TRAIT_MODE } from '../subworkflows/RERCONVERGE/rer_trait_mode'
include { RER_CONT } from '../subworkflows/RERCONVERGE/rer_cont'
include { RER_BIN } from '../subworkflows/RERCONVERGE/rer_bin'
include { RER_TABLES } from '../subworkflows/RERCONVERGE/rer_tables'
include { RER_REPORT } from '../subworkflows/RERCONVERGE/rer_report'
include { RER_GENE_LISTS } from '../subworkflows/RERCONVERGE/rer_gene_lists'
include { rerSingle; rerRequired; rerStage } from '../subworkflows/RERCONVERGE/rer_contract'
include { RER_FCS_REPORT } from '../subworkflows/ENRICHMENT/fcs.nf'
include { FCS_RER as FCS_COMPUTE } from '../subworkflows/ENRICHMENT/fcs_inputs.nf'

workflow RER_MAIN {
    take:
    traitfile_input
    pp_top_ch       // retained public signature; unused by RER
    pp_bottom_ch
    taxid_map_input

    main:
    def toolString = params.rer_tool ? params.rer_tool.toString().trim() : ''
    def steps = toolString && toolString != 'false' ? toolString.split(',').collect { it.trim() } : []
    assert steps.every { it in ['build_trait','build_tree','build_matrix','continuous','binary'] } : 'Unknown RER step'
    assert steps.unique(false).size()==steps.size() : 'Duplicate RER steps'
    assert (params.rer_trait_mode ?: 'auto') in ['auto','binary','continuous'] : 'Invalid rer_trait_mode'
    assert params.traitname && !(params.traitname.toString() =~ /[\/\t\r\n]/).find() : 'Invalid RER traitname'
    def correlation = steps.contains('continuous') || steps.contains('binary')
    assert !(params.rer_continuous_file && params.rer_binary_file) : 'Choose exactly one precomputed RER result'
    assert !(correlation && (params.rer_continuous_file || params.rer_binary_file)) : 'Cannot combine new and precomputed RER results'
    assert !(correlation && params.rer_perms_file) : 'New RER inference produces its own null; remove rer_perms_file'

    def raw = Channel.empty()
    if (steps.contains('build_trait') || steps.contains('build_tree')) {
        raw = rerSingle((traitfile_input ?: Channel.empty()).ifEmpty {
            assert params.my_traits : 'RER requires my_traits'
            file(params.my_traits)
        }, 'trait table')
    }
    def trait = Channel.empty()
    if (steps.contains('build_trait')) trait = RER_TRAIT(raw,rerStage('trait')).polished
    else if (steps.contains('build_matrix') || correlation) trait = rerRequired(params.trait_out,'trait_out')

    def trees = Channel.empty()
    if (steps.contains('build_tree')) {
        def taxid = (taxid_map_input ?: Channel.empty()).toList().map { xs ->
            assert xs.size() <= 1 : 'RER: expected at most one TaxID table'
            def f = xs ? xs[0] : file('NO_FILE')
            assert f.name=='NO_FILE' || f.isFile() : 'RER: TaxID table not found'
            f
        }
        trees = RER_TREES(raw,rerRequired(params.gene_trees,'gene_trees'),taxid,rerStage('trees')).master_tree
    } else if (steps.contains('build_matrix') || correlation) trees = rerRequired(params.trees_out,'trees_out')

    def matrix = Channel.empty()
    if (steps.contains('build_matrix')) matrix = RER_MATRIX(trait,trees,rerStage('matrix')).matrix
    else if (correlation) matrix = rerRequired(params.matrix_out,'matrix_out')

    def results = Channel.empty()
    def nulls = Channel.empty()
    if (correlation) {
        def mode = RER_TRAIT_MODE(trait,rerStage('mode')).mode.map { it.text.trim() }
        def routed = trait.combine(mode).branch {
            continuous: it[1]=='continuous'
            binary: it[1]=='binary'
        }
        def cont = RER_CONT(routed.continuous.map { f,t -> f },trees,matrix,rerStage('continuous'))
        def bin = RER_BIN(routed.binary.map { f,t -> f },trees,matrix,rerStage('binary'))
        results = cont.continuous_output.mix(bin.binary_output)
        nulls = cont.perms_output.mix(bin.perms_output)
    } else if (params.rer_continuous_file || params.rer_binary_file) {
        results = rerRequired(params.rer_continuous_file ?: params.rer_binary_file,'precomputed result')
    }
    if (!correlation && params.rer_perms_file) nulls = rerRequired(params.rer_perms_file,'rer_perms_file')
    def null_input = nulls.ifEmpty(file('NO_FILE'))
    def tables = RER_TABLES(results,null_input,rerStage('tables'))
    RER_REPORT(tables.validated,rerStage('report'))
    def gene_scores = params.rer_gene_scores ? rerRequired(params.rer_gene_scores,'rer_gene_scores') : Channel.value(file('NO_FILE'))
    def lists = RER_GENE_LISTS(tables.summary_tsv,gene_scores,rerStage('lists'))
    def background = lists.gene_lists.flatten().filter { it.name=='background.txt' }.collect()
    def interest = lists.gene_lists.flatten().filter { it.name!='background.txt' }.collect()
    // Release the null only after its association, roster and checksums were validated.
    def verified_null = nulls.combine(tables.summary_tsv).map { p,s -> p }
    if (params.enrichment && !params.scoring) {
        def fcs_null = verified_null.ifEmpty(file('NO_FILE'))
        def fcs = FCS_COMPUTE(lists.fcs_stats,background,fcs_null)
        RER_FCS_REPORT(Channel.value('rerconverge'),lists.fcs_stats,background,
            Channel.value("11.FCS_rer_${params.traitname}"),fcs_null,
            Channel.value(file('NO_FILE')),fcs.enrich_file)
    }
    emit:
    summary_tsv = tables.summary_tsv
    perms = verified_null
    gene_lists_bg = background
    gene_lists_interest = interest
}
