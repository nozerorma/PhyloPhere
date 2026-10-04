// Runs the real CAAS_CORE (batching + CAAS_CORE_BATCHED) on the given alignments or perm-discovery exports and
// writes the work directories of the shard directories to <outdir>/pos_detail_dirs.txt and those of the b_0
// observed directories to <outdir>/b0_observed_dirs.txt; used as the main.nf of a project that links the tree
// under test (see test_wiring.py).
include { CAAS_CORE } from './subworkflows/CT/caas_permulation'

workflow {
    def alignments = params.mini_alignments
        ? Channel.fromList(params.mini_alignments.split(',').collect { f -> tuple(file(f).baseName, file(f)) })
        : Channel.empty()
    def reuse = params.mini_reuse
        ? Channel.value(params.mini_reuse.split(',').collect { f -> file(f) })
        : Channel.empty()
    CAAS_CORE(alignments, reuse, Channel.value(file(params.mini_cfg)), Channel.value(file(params.mini_resample)),
              Channel.value(file(params.mini_tree)), Channel.value(file(params.mini_fop_pairs ?: 'NO_FOP_PAIRS')), Channel.value(file('NO_FILE')))
    CAAS_CORE.out.pos_detail.map { d -> d.toString() }.collectFile(name: 'pos_detail_dirs.txt', storeDir: params.outdir, newLine: true)
    CAAS_CORE.out.b0_observed.map { d -> d.toString() }.collectFile(name: 'b0_observed_dirs.txt', storeDir: params.outdir, newLine: true)
}
