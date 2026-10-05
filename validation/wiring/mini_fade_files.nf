// Prints the alignment tuples and the toy genes of the FADE selection functions, next to the sample CT draws from the same
// directory; used as the main.nf of a project that links the tree under test (see test_wiring.py).
include { ali_tuples_from_dir; fade_toy_genes } from './subworkflows/SELECTION/selection_prep'
include { listAlignmentFiles; sampleAlignmentFiles } from './subworkflows/CT/ct_alignment_files'

workflow {
    def n = params.mini_n as int
    def seed = params.mini_seed as long
    def out = [
        listed : ali_tuples_from_dir(params.mini_dir, null).collect { t -> t[1].name },
        toy    : fade_toy_genes(params.mini_dir, n, seed),
        ct_toy : sampleAlignmentFiles(listAlignmentFiles(params.mini_dir), n, seed).collect { f -> f.name.tokenize('.')[0] },
        wanted : ali_tuples_from_dir(params.mini_dir, ['GENE03', 'GENE01'] as Set).collect { t -> t[1].name },
    ]
    file(params.out).text = groovy.json.JsonOutput.toJson(out)
}
