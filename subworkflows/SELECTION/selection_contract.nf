#!/usr/bin/env nextflow
include { fingerprintRoots } from '../common/content_hash.nf'

def selectionQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def selectionHash(sources) {
    sources.each { source -> assert source.isFile() : "SELECTION_PREP: file not found: ${source}" }
    return fingerprintRoots(sources, true)
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
