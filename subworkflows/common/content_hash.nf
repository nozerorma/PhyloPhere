#!/usr/bin/env nextflow

// Source fingerprints and data metadata for task cache keys. One entry hashes as: name bytes, a 0 byte, the content,
// then a 0 byte unless nulAfterContent is false. A directory contributes its files recursively, sorted
// by path and named relative to the directory; a file contributes itself under its own name.
// Policy change intentionally invalidates older deep-hash task keys once.
// Callers assert that their inputs exist; a missing root contributes nothing.

// Hash small implementation sources; identify data by the same metadata policy
// as ordinary Nextflow caching. Do not read multi-GB resources at workflow startup.
def fingerprintFile(md, java.nio.file.Path f) {
    def source = f.toString() ==~ /.*\/(?:subworkflows|workflows|bin)\/.*\.(?:py|R|r|nf|Rmd|qmd|sh|bf)/
    if (source) {
        f.toFile().withInputStream { stream ->
            new java.security.DigestInputStream(stream, md).transferTo(java.io.OutputStream.nullOutputStream())
        }
    } else {
        def stat = java.nio.file.Files.readAttributes(f, java.nio.file.attribute.BasicFileAttributes)
        md.update((stat.size().toString()+':'+stat.lastModifiedTime().toString()).getBytes('UTF-8'))
    }
}

def fingerprintEntry(md, String name, java.nio.file.Path f, boolean nulAfterContent) {
    md.update(name.getBytes('UTF-8'))
    md.update(0 as byte)
    fingerprintFile(md, f)
    if (nulAfterContent) md.update(0 as byte)
}

// Hash one root (file or directory) into md.
def fingerprintRoot(md, java.nio.file.Path root, boolean nulAfterContent) {
    if (!root.toFile().exists()) return
    if (!root.toFile().isDirectory()) {
        fingerprintEntry(md, root.name, root, nulAfterContent)
        return
    }
    def files = []
    root.toFile().eachFileRecurse(groovy.io.FileType.FILES) { files << it.toPath() }
    files.sort { it.toString() }.each { f ->
        fingerprintEntry(md, root.relativize(f).toString(), f, nulAfterContent)
    }
}

// SHA-256 hex digest of the roots, in the order given (the caller sorts when its layout needs it).
def fingerprintRoots(List roots, boolean nulAfterContent) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    roots.each { fingerprintRoot(md, it, nulAfterContent) }
    md.digest().encodeHex().toString()
}
