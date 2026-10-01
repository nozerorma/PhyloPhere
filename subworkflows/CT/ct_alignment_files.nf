#!/usr/bin/env nextflow

/*
 * Alignment files of a run, shared by CT and the standalone permulation null.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

// Alignment files of a directory in file-name order, without the annotation tables and logs that can sit
// beside them. An absent directory gives an empty list.
def listAlignmentFiles(dirPath) {
    def dir = file(dirPath as String)
    if (!dir.exists()) return []
    def files = dir.listFiles()?.findAll { f -> f.isFile() && !f.name.matches('.*\\.txt$|.*\\.tsv$|.*\\.csv$|.*\\.log$|.*\\.map$') } ?: []
    return files.sort(false) { f -> f.name }
}

// toy_mode: n alignments chosen by a seeded shuffle of the name-sorted list, so the same seed picks the same
// alignments whatever order the filesystem lists them in. Batch composition (and with it -resume cache hits)
// follows the order of the list.
def sampleAlignmentFiles(files, n, seed) {
    def shuffled = new ArrayList(files.sort(false) { f -> f.name })
    Collections.shuffle(shuffled, new Random(seed as long))
    return shuffled.take(n as int)
}
