// Prints the resolved defaults under test; used as the main.nf of a project that links the tree under test (see test_wiring.py).
workflow {
    println "PARAMS filter_maxcaas=${params.filter_maxcaas} caas_map_dir=[${params.caas_map_dir}]"
}
