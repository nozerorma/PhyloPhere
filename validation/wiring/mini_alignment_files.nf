// Prints the alignment list and the toy sample of the functions CT and the standalone null share, for the input
// orders given in params.mini_orders; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { listAlignmentFiles; sampleAlignmentFiles } from './subworkflows/CT/ct_alignment_files'

workflow {
    def listed = listAlignmentFiles(params.mini_dir)
    def out = [listed: listed.collect { f -> f.name }]
    params.mini_orders.split(';').each { order ->
        def names = order.split(',') as List
        def files = names.collect { n -> listed.find { f -> f.name == n } }
        out[order] = sampleAlignmentFiles(files, params.mini_n as int, params.mini_seed as long).collect { f -> f.name }
    }
    file(params.out).text = groovy.json.JsonOutput.toJson(out)
}
