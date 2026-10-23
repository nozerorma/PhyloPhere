// Directed integration on a saved WP5 run, without replaying the permulation core.
// Copy to the snapshot root and rewrite ../ imports to ./ (prepare_wp7_snapshot.sh).
nextflow.enable.dsl=2
include { CAAS_CORE_MERGE_REDUCE; CAAS_CORE_MERGE_REDUCE as CAAS_CORE_MERGE_REDUCE_STRICT } from '../subworkflows/CT/caas_permulation.nf'
include { SCORING } from '../workflows/scoring.nf'
include { FCS_CAAS } from '../subworkflows/ENRICHMENT/fcs_inputs.nf'
include { SCORING_FCS_REPORT } from '../subworkflows/ENRICHMENT/fcs.nf'
workflow {
    def saved = params.wp7_saved
    def universe = Channel.value(file("${saved}/universes/cleaned_background_main_marine.txt"))
    def blocks = Channel.fromPath("${saved}/caas_permulation/blocks/blk_*",type:'dir',checkIfExists:true).collect()
    def reduced = CAAS_CORE_MERGE_REDUCE(blocks, universe, Channel.value(file("${saved}/caas_permulation/resample_perms.tab")))
    CAAS_CORE_MERGE_REDUCE_STRICT(Channel.fromPath("${saved}/strict/caas_permulation/blocks/blk_*",type:'dir',checkIfExists:true).collect(),
       universe, Channel.value(file("${saved}/caas_permulation/resample_perms.tab")))
    def validated = blocks.map { bs -> [blocks:bs] }.combine(reduced.metadata).map { wrapper,meta -> wrapper.blocks }
    SCORING(
        Channel.value(file("${saved}/postproc/gene_filtering/filtered_discovery.tsv")),
        null, null,
        Channel.value(file("${saved}/rerconverge/rer_results/rerconverge_summary_marine.tsv")),
        Channel.fromPath("${saved}/accumulation/**/*.csv"),
        null, null, null, universe, null,
        reduced.perms, reduced.null_positions, validated, reduced.counts,
        reduced.cycle_summary.combine(reduced.cycle_hist), reduced.pos_sample, reduced.pos_quantiles,
        null,
        Channel.value(file("${saved}/caastools/permulation_manifest.tsv")),
        Channel.value(file("${saved}/caastools/permulation_harvest.tsv")),
        Channel.value(file("${saved}/caastools/background.output")), reduced.cycle_roster
    )
    FCS_CAAS(SCORING.out.fcs_stats,universe,SCORING.out.caas_perms,validated)
    SCORING_FCS_REPORT(SCORING.out.fcs_stats,universe,SCORING.out.gene_lists,FCS_CAAS.out.enrich_file)
    FCS_CAAS.out.enrich_file.subscribe { f -> println "WP7_FCS=${f}" }
}
