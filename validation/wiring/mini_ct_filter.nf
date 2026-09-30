// Runs the real CT_FILTER process on one input; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { CT_FILTER } from './subworkflows/CT_POSTPROC/ctpp_clustfilter'

workflow {
    CT_FILTER(Channel.of(tuple('filter', 3, 0.7d, file(params.mini_input))))
}
