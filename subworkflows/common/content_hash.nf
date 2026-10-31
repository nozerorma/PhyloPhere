#!/usr/bin/env nextflow

// Content fingerprints for task cache keys. One entry hashes as: name bytes, a 0 byte, the content,
// then a 0 byte unless nulAfterContent is false. A directory contributes its files recursively, sorted
// by path and named relative to the directory; a file contributes itself under its own name.
// The byte layout is fixed: changing it changes every cache key that depends on it.
// Callers assert that their inputs exist; a missing root contributes nothing.

// Stream one file through the digest in chunks, so large files are never loaded whole.
def fingerprintFile(md, java.nio.file.Path f) {
    f.toFile().withInputStream { stream ->
        new java.security.DigestInputStream(stream, md).transferTo(java.io.OutputStream.nullOutputStream())
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
