// Runs the real CAAS_CORE_MERGE on the given shard directories (or one legacy file) and writes the work directory of
// the merged shard directory to <outdir>/pos_detail_dir.txt when there is one; used as the main.nf of a project that
// links the tree under test (see test_wiring.py).
include { CAAS_CORE_MERGE } from './subworkflows/CT/caas_permulation'

workflow {
    def b0 = params.mini_b0 ? params.mini_b0.split(',').collect { f -> file(f) } : file('NO_B0_OBSERVED')
    def merged = CAAS_CORE_MERGE(Channel.value(params.mini_details.split(',').collect { f -> file(f) }),
                                 Channel.value(file('NO_FILE')), Channel.value(file(params.mini_lengths)),
                                 Channel.value(b0), Channel.value(file(params.mini_design ?: 'NO_DESIGN')))
    merged.pos_detail.map { d -> d.toString() }.collectFile(name: 'pos_detail_dir.txt', storeDir: params.outdir, newLine: true)
    merged.master.mix(merged.discovery, merged.background, merged.background_genes, merged.meta_caas).map { f -> f.toString() }
        .collectFile(name: 'contract_paths.txt', storeDir: params.outdir, newLine: true)
}
