#!/usr/bin/env nextflow
// Directed regeneration: published seeds/backgrounds -> fresh DOMINO -> AMI.
// Place this entry point at the repository root; it intentionally does not import main.nf.
nextflow.enable.dsl = 2

include { DOMINO_MODULES as REGEN_CAAS } from './subworkflows/ENRICHMENT/domino.nf'
include { DOMINO_MODULES as REGEN_FADE } from './subworkflows/ENRICHMENT/domino.nf'
include { DOMINO_MODULES as REGEN_RER } from './subworkflows/ENRICHMENT/domino.nf'
include { SCORING_AMI_REPORT } from './subworkflows/ENRICHMENT/scoring_enrichment.nf'

params.ami_source = null

workflow {
    if (!params.ami_source) error 'Provide --ami_source pointing to the published run'
    if (!params.traitname) error 'Provide --traitname for the published run'
    def source = file(params.ami_source, checkIfExists: true)
    if (!source.isDirectory()) error "Not a published results directory: ${source}"
    if (source.toAbsolutePath().normalize() == file(params.outdir).toAbsolutePath().normalize())
        error 'Use a separate --outdir to review regenerated AMI before publication'

    def required = { String rel -> file("${source}/${rel}", checkIfExists: true) }
    def trait = params.traitname.toString()
    def seeds = required('scoring/gene_lists')
    ['global','top','bottom'].each { required("scoring/gene_lists/active_${it}.tsv") }
    def scores = required('scoring/gene_scores.tsv')
    def metadata = required('scoring/gene_evidence_metadata.json')
    def stamp = new groovy.json.JsonSlurper().parseText(metadata.text)
    if (stamp.gene_stat != 'max') error 'Published gene evidence must record gene_stat=max'
    def caas_bg = required("universes/cleaned_background_main_${trait}.txt")
    def fade_bg = required("universes/fade_background_${trait}.txt")
    def rer_bg = required("universes/rer_background_${trait}.txt")
    def fade_seeds = [required('selection/fade/top/gene_lists/fade_top_significant.txt'),
                      required('selection/fade/bottom/gene_lists/fade_bottom_significant.txt')]
    def rer_seeds = ['rer_significant','rer_accelerating','rer_decelerating'].collect {
        required("rerconverge/gene_lists/${it}.txt")
    }
    // The shared list preparer builds FADE Global as the union of these two sides.
    // Each method retains its own published background and independently rebuilt network.
    REGEN_CAAS(Channel.value(caas_bg), Channel.value(seeds),
        params.domino_network_score_thr, params.domino_slice_thr,
        params.domino_module_thr, params.string_db_dir ?: '', 'caas')
    REGEN_FADE(Channel.value(fade_bg), Channel.value(fade_seeds),
        params.domino_network_score_thr, params.domino_slice_thr,
        params.domino_module_thr, params.string_db_dir ?: '', 'fade')
    REGEN_RER(Channel.value(rer_bg), Channel.value(rer_seeds),
        params.domino_network_score_thr, params.domino_slice_thr,
        params.domino_module_thr, params.string_db_dir ?: '', 'rer')

    SCORING_AMI_REPORT(
        Channel.value(seeds), Channel.value(caas_bg), Channel.value(scores), Channel.value(metadata),
        REGEN_CAAS.out.network_sif, REGEN_CAAS.out.modules_dir, REGEN_CAAS.out.edge_scores,
        Channel.value(fade_seeds), Channel.value(fade_bg),
        REGEN_FADE.out.network_sif, REGEN_FADE.out.modules_dir, REGEN_FADE.out.edge_scores,
        Channel.value(rer_seeds), Channel.value(rer_bg),
        REGEN_RER.out.network_sif, REGEN_RER.out.modules_dir, REGEN_RER.out.edge_scores,
        true, true
    )
}
