// Runs the real CAAS_CORE_MERGE on the given shard directories (or one legacy file) and writes the work directory of
// the merged shard directory to <outdir>/pos_detail_dir.txt when there is one; used as the main.nf of a project that
// links the tree under test (see test_wiring.py).
include { CAAS_CORE_MERGE } from './subworkflows/CT/caas_permulation'

workflow {
    def merged = CAAS_CORE_MERGE(Channel.value(params.mini_details.split(',').collect { f -> file(f) }),
                                 Channel.value(file('NO_FILE')), Channel.value(file(params.mini_lengths)),
                                 Channel.value(file(params.mini_labelings ?: 'NO_LABELINGS')))
    merged.pos_detail.map { d -> d.toString() }.collectFile(name: 'pos_detail_dir.txt', storeDir: params.outdir, newLine: true)
}
