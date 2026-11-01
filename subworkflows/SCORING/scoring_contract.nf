// Content fingerprints are interpolated into each task script before cache lookup.
def scoringHash(paths) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    md.update(('integrity=' + (System.getenv('PHYLOPHERE_VERIFY_CHECKSUMS') ?: '0')).getBytes('UTF-8'))
    paths.sort { it.toString() }.each { p ->
        assert p.isFile() : "Missing SCORING source: ${p}"
        md.update(p.name.getBytes('UTF-8'))
        md.update('\u0000'.getBytes('UTF-8'))
        md.update(p.bytes)   // sources are small scripts, read whole
    }
    md.digest().encodeHex().toString()
}
def scoringQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
