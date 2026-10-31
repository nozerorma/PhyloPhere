// VEP content fingerprints are explicit inputs, including copied Python sources.
def vepQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def vepHash(paths) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    paths.collect { file(it) }.sort { it.toString() }.each { root ->
        assert root.exists() : "Missing VEP input: ${root}"
        def files = []
        if (root.isDirectory()) {
            root.toFile().eachFileRecurse(groovy.io.FileType.FILES) { f -> files << f.toPath() }
        } else { files << root }
        files.sort { it.toString() }.each { p ->
            def name = root.isDirectory() ? root.relativize(p).toString() : root.name
            md.update(name.getBytes('UTF-8')); md.update((byte)0)
            p.toFile().withInputStream { stream ->
                byte[] buffer = new byte[1048576]
                int n
                while ((n = stream.read(buffer)) != -1) md.update(buffer,0,n)
            }
        }
    }
    md.digest().encodeHex().toString()
}

def vepSources(names) {
    vepHash(['vep_contract.nf'] .collect { "${baseDir}/subworkflows/VEP/${it}" } +
        names.collect { "${baseDir}/subworkflows/VEP/${it}" })
}
