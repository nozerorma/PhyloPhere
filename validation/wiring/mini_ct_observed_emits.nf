// Runs the real CT_OBSERVED on the given files and writes the design and tree it emits to <outdir>/emits.txt;
// used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CT_OBSERVED } from './workflows/ct_observed'

workflow {
    def obs = CT_OBSERVED(Channel.value(file(params.mini_discovery)), Channel.value(file(params.mini_design)),
                          Channel.value(file(params.mini_tree)), null)
    obs.design.map { f -> "design ${f}".toString() }.mix(obs.tree.map { f -> "tree ${f}".toString() })
        .collectFile(name: 'emits.txt', storeDir: params.outdir, newLine: true, sort: true)
}
