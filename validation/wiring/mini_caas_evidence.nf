// Runs the real CAAS_EVIDENCE on the given files and writes the paths of the files it emits to <outdir>/evidence_paths.txt;
// used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CAAS_EVIDENCE } from './subworkflows/CT_DISAMBIGUATION/ct_evidence'

workflow {
    def ev = CAAS_EVIDENCE(Channel.value(file(params.mini_discovery)), Channel.value(file(params.mini_scores)),
                           Channel.value(file(params.mini_design)), Channel.value(file(params.mini_tree)))
    ev.evidence_tsv.mix(ev.top_positions).map { f -> f.toString() }
        .collectFile(name: 'evidence_paths.txt', storeDir: params.outdir, newLine: true)
}
