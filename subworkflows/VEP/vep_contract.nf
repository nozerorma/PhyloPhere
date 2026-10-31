// VEP content fingerprints are explicit inputs, including copied Python sources.
include { fingerprintRoots } from '../common/content_hash.nf'
def vepQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def vepHash(paths) {
    // VEP entries have no separator after the content (nulAfterContent = false); this is the existing cache layout.
    def roots = paths.collect { file(it) }.sort { it.toString() }
    roots.each { root -> assert root.exists() : "Missing VEP input: ${root}" }
    return fingerprintRoots(roots, false)
}

def vepSources(names) {
    vepHash(['vep_contract.nf'] .collect { "${baseDir}/subworkflows/VEP/${it}" } +
        names.collect { "${baseDir}/subworkflows/VEP/${it}" })
}
