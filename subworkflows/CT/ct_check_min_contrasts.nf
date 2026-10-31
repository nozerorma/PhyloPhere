#!/usr/bin/env nextflow
// ct_check_min_contrasts.nf — Gate that stops a trait with too few foreground contrasts before CT.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CHECK_MIN_CONTRASTS: counts the foreground species (column 2 == 1) of the
 *  caastools traitfile produced by CONTRAST_ALGORITHM and compares the count with
 *  params.min_contrasts (3 when unset).
 *
 *  Threshold not met: a sentinel, low_contrasts.skip (trait, count, minimum), is
 *  published to params.outdir and no traitfile is emitted, so the processes that
 *  consume them do not run. main.nf stops the run gracefully (exit 0) when it sees
 *  the sentinel, and the single-phenotype launch scripts check for the file to move on
 *  to the next phenotype.
 *
 *  Threshold met: the traitfile, the permulation traitfile and the traitfile
 *  directory are copied through under *_ok names and no sentinel is written.
 *
 *  Consumes:  traitfile, permulation traitfile, traitfile directory
 *             (CONTRAST_ALGORITHM, through workflows/contrast_selection.nf)
 *  Produces:  traitfile_ok.tab, permulation_traitfile_ok.tab, traitfiles_ok_dir/
 *             or low_contrasts.skip
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Minimum-contrast gate ────────────────────────────────────────────────────

process CHECK_MIN_CONTRASTS {
    tag "CHECK_MIN_CONTRASTS"
    label 'process_discovery'   // lightest available label; trivially fast

    input:
    path traitfile
    path permulation_traitfile
    path trait_dir, stageAs: 'trait_dir_in'

    output:
    path "gate_status.json", emit: status
    path "traitfile_ok.tab",             emit: traitfile_out,             optional: true
    path "permulation_traitfile_ok.tab", emit: permulation_traitfile_out, optional: true
    path "traitfiles_ok_dir",            emit: trait_dir_out,             optional: true
    path "low_contrasts.skip",           emit: skip_flag,                 optional: true

    script:
    def arguments = groovy.json.JsonOutput.toJson([
        traitfile: traitfile.toString(), permulation_traitfile: permulation_traitfile.toString(),
        trait_dir: trait_dir.toString(), min_contrasts: params.min_contrasts != null ? params.min_contrasts : 3,
        traitname: params.traitname ?: 'unknown', multi_hypothesis: params.multi_hypothesis == true
    ])
    """
    cp "${baseDir}/subworkflows/CT/local/scripts/check_contrasts.py" .
    cat > gate_inputs.json <<'PHYLOPHERE_GATE_INPUTS'
${arguments}
PHYLOPHERE_GATE_INPUTS
    python check_contrasts.py gate_inputs.json
    """
}

// Reconcile on the Nextflow head, including cache hits. Optional publishDir
// outputs cannot remove an old sentinel, and publication can race main.nf's exit.
// Every public output depends on this reconciliation completing first.
workflow GATE_CONTRASTS {
    take:
    traitfile
    permulation_traitfile
    trait_dir

    main:
    result = CHECK_MIN_CONTRASTS(traitfile, permulation_traitfile, trait_dir)
    reconciled = result.status.map { status_file ->
        def status = new groovy.json.JsonSlurper().parse(status_file.toFile())
        def destination = file(params.outdir).toAbsolutePath()
        java.nio.file.Files.createDirectories(destination)
        def sentinel = destination.resolve('low_contrasts.skip')
        if (status.skipped) {
            def temporary = java.nio.file.Files.createTempFile(destination, '.low_contrasts-', '.tmp')
            try {
                java.nio.file.Files.copy(status_file.parent.resolve('low_contrasts.skip'), temporary,
                    java.nio.file.StandardCopyOption.REPLACE_EXISTING)
                java.nio.file.Files.move(temporary, sentinel, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            } finally {
                java.nio.file.Files.deleteIfExists(temporary)
            }
        } else {
            java.nio.file.Files.deleteIfExists(sentinel)
        }
        true
    }

    emit:
    traitfile_out = result.traitfile_out.combine(reconciled).map { path, ready -> path }
    permulation_traitfile_out = result.permulation_traitfile_out.combine(reconciled).map { path, ready -> path }
    trait_dir_out = result.trait_dir_out.combine(reconciled).map { path, ready -> path }
    skip_flag = result.skip_flag.combine(reconciled).map { path, ready -> path }
}
