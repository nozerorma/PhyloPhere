// Runs the real CAAS_PERMS_PREP (SUBSET_RESAMPLE_PERMS) and writes the work paths of the resample subset and of
// the FOP pairs to <outdir>/subset_paths.txt; used as the main.nf of a project that links the tree under test
// (see test_wiring.py).
include { CAAS_PERMS_PREP } from './subworkflows/CT/caas_permulation'

workflow {
    def prep = CAAS_PERMS_PREP(Channel.value(file(params.mini_cfg)), Channel.value(file(params.mini_resample)))
    prep.resample_subset.mix(prep.fop_pairs).map { f -> f.toString() }
        .collectFile(name: 'subset_paths.txt', storeDir: params.outdir, newLine: true)
}
