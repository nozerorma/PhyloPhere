// Calls the post-processing helper functions of the real modules; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { clusterParameterGrid; clusterFileSuffix } from './subworkflows/CT_POSTPROC/ctpp_clustfilter'
include { caasPostprocArgs; caasPostprocOn } from './subworkflows/CT/caas_permulation'

// Doubles, as in the workflow (params.filter_maxcaas.toDouble()): a Groovy 0.29 literal is a BigDecimal and would give 29.
workflow {
    println "GRID_IN=" + clusterParameterGrid([2, 3, 4, 10], [0.6, 0.7, 0.8], 3, 0.7)
    println "GRID_OUT=" + clusterParameterGrid([2, 3, 4], [0.6, 0.7], 5, 0.65)
    println "SUFFIX_70=" + clusterFileSuffix(3, 0.7d)
    println "SUFFIX_29=" + clusterFileSuffix(2, 0.29d)
    println "ARGS=[" + caasPostprocArgs(file('genes.tsv')) + "]"
    println "ARGS_SENTINEL=[" + caasPostprocArgs(file('NO_FILE')) + "]"
    println "ON=" + caasPostprocOn()
}
