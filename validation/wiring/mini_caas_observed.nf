// Runs the real CAAS_OBSERVED (a discovery.tab scored without a replay) and writes the paths of the files it emits
// to <outdir>/observed_paths.txt; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CAAS_OBSERVED } from './subworkflows/CT_DISAMBIGUATION/ct_observed'

workflow {
    def obs = CAAS_OBSERVED(Channel.value(file(params.mini_discovery)), Channel.value(file(params.mini_design)),
                            Channel.value(file(params.mini_tree)), Channel.value(file(params.mini_hyp_pairs)))
    obs.master_csv.mix(obs.global_meta_caas).map { f -> f.toString() }
        .collectFile(name: 'observed_paths.txt', storeDir: params.outdir, newLine: true)
}
