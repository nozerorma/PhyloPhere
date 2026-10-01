// Runs the real CAAS_PERMS_PREP (subset + batched perm-replay) on the given alignments; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CAAS_PERMS_PREP } from './subworkflows/CT/caas_permulation'

workflow {
    def alignments = Channel.fromList(params.mini_alignments.split(',').collect { f -> tuple(file(f).baseName, file(f)) })
    CAAS_PERMS_PREP(alignments, file(params.mini_cfg), file(params.mini_resample))
}
