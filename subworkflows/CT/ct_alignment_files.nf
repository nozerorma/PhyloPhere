#!/usr/bin/env nextflow
// ct_alignment_files.nf — List and subsample the alignment files of a run.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  Groovy helpers (no process) shared by the CT workflow, the standalone permulation
 *  null in main.nf and the SELECTION preparation.
 *
 *  Consumes:  an alignment directory (params.alignment)
 *  Produces:  a name-sorted list of alignment files, optionally a seeded subset of it
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Helpers ──────────────────────────────────────────────────────────────────

// Alignment files of a directory in file-name order, without the annotation tables and logs that can sit
// beside them (.txt, .tsv, .csv, .log, .map). An absent directory gives an empty list.
def listAlignmentFiles(dirPath) {
    def dir = file(dirPath as String)
    if (!dir.exists()) return []
    def files = dir.listFiles()?.findAll { f -> f.isFile() && !f.name.matches('.*\\.txt$|.*\\.tsv$|.*\\.csv$|.*\\.log$|.*\\.map$') } ?: []
    return files.sort(false) { f -> f.name }
}

// Used by toy_mode: n alignments chosen by a seeded shuffle of the name-sorted list, so the same seed picks the
// same alignments whatever order the filesystem lists them in. Batch composition, and with it the -resume cache
// hits, follows the order of the returned list.
def sampleAlignmentFiles(files, n, seed) {
    def shuffled = new ArrayList(files.sort(false) { f -> f.name })
    Collections.shuffle(shuffled, new Random(seed as long))
    return shuffled.take(n as int)
}
