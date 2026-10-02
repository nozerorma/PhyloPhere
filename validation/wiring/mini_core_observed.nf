// Runs the real CAAS_CORE_OBSERVED on the given b_0 directories and writes the paths of the files it emits to
// <outdir>/contract_paths.txt; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CAAS_CORE_OBSERVED } from './subworkflows/CT/caas_permulation'

workflow {
    def obs = CAAS_CORE_OBSERVED(Channel.value(params.mini_b0.split(',').collect { f -> file(f) }), Channel.value(file(params.mini_design)))
    obs.discovery.mix(obs.background, obs.background_genes, obs.master_csv, obs.global_meta_caas)
        .map { f -> f.toString() }.collectFile(name: 'contract_paths.txt', storeDir: params.outdir, newLine: true)
}
