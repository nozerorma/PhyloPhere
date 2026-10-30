nextflow.enable.dsl=2
include { SCORING_EVIDENCE } from './subworkflows/SCORING/scoring_evidence.nf'
include { FCS_COMPUTE_CALIBRATED } from './subworkflows/ENRICHMENT/fcs_blocks.nf'
include { SCORING } from './workflows/scoring.nf'
include { ENRICHMENT } from './workflows/enrichment.nf'
workflow {
    def input = params.fixture
    def universe = Channel.value(file("${input}/universe.txt"))
    SCORING_EVIDENCE(Channel.value(file("${input}/genes.tsv")),Channel.value(file("${input}/positions.tsv")),
        Channel.value(file("${input}/stats.tsv")),universe,Channel.value([]),Channel.value(file("${input}/null.rds")))
    SCORING_EVIDENCE.out.gene_evidence.subscribe { println "EVIDENCE=${it}" }
    if(params.fcs_enabled) {
        FCS_COMPUTE_CALIBRATED(universe, SCORING_EVIDENCE.out.calibration)
        FCS_COMPUTE_CALIBRATED.out.enrich_file.subscribe { println "CALIBRATED_FCS=${it}" }
    }
}
