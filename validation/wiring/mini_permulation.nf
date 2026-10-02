// Runs the real CAAS_PERMULATION workflow (batches + merge) on the given alignments and writes the paths of the
// observed contract files it emits to <outdir>/contract_paths.txt; used as the main.nf of a project that links the
// tree under test (see test_wiring.py).
include { CAAS_PERMULATION } from './subworkflows/CT/caas_permulation'

workflow {
    def alignments = Channel.fromList(params.mini_alignments.split(',').collect { f -> tuple(file(f).baseName, file(f)) })
    def perm = CAAS_PERMULATION(alignments, Channel.empty(), Channel.value(file(params.mini_cfg)), Channel.value(file(params.mini_resample)),
                                Channel.value(file(params.mini_tree)), Channel.value(file('NO_FILE')), Channel.value(file(params.mini_fop_pairs)),
                                Channel.value(file(params.mini_lengths)), Channel.value('NO_GATE'))
    perm.discovery.mix(perm.background, perm.background_genes, perm.master, perm.global_meta_caas)
        .map { f -> f.toString() }.collectFile(name: 'contract_paths.txt', storeDir: params.outdir, newLine: true)
}
