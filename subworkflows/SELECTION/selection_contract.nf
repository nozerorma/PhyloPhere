#!/usr/bin/env nextflow

def selectionQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def selectionHash(sources) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    sources.each { source ->
        assert source.isFile() : "SELECTION_PREP: file not found: ${source}"
        md.update(source.name.getBytes('UTF-8'))
        md.update((byte)0)
        source.toFile().withInputStream { stream ->
            byte[] buffer = new byte[1048576]
            int n
            while ((n = stream.read(buffer)) != -1) md.update(buffer, 0, n)
        }
        md.update((byte)0)
    }
    md.digest().encodeHex().toString()
}

def selectionSingle(files, label) {
    assert files.size() == 1 : "SELECTION_PREP: expected one ${label}, got ${files.size()}"
    assert files[0].isFile() : "SELECTION_PREP: ${label} not found: ${files[0]}"
    files[0]
}

def selectionPositive(value, label) {
    assert value.toString() ==~ /[1-9][0-9]*/ : "SELECTION_PREP: ${label} must be a positive integer"
    value as int
}
