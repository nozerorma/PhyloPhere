#!/usr/bin/env nextflow

/*
#
# PHYLOPHERE — Minimum-contrast gate
#
# Checks that the caastools traitfile produced by CONTRAST_ALGORITHM contains
# at least params.min_contrasts (default 3) foreground species (column 2 == 1).
#
# When the threshold is NOT met:
#   • A plain-text sentinel  low_contrasts.skip  is published to ${params.outdir}.
#   • Neither traitfile nor boot_traitfile is emitted downstream.
#   • All subsequent CT / signification / disambiguation / … processes are
#     silently skipped because their input channels never receive a value.
#   • The Nextflow run exits 0 (no error).
#   • The bash orchestrators (run_phenotypes.sh, test_stress.sh) detect the
#     sentinel and continue to the next phenotype.
#
# When the threshold IS met:
#   • Both traitfiles are emitted unchanged (as traitfile_ok.tab /
#     boot_traitfile_ok.tab) and downstream processing proceeds normally.
#   • No sentinel file is written.
*/

process CHECK_MIN_CONTRASTS {
    tag "CHECK_MIN_CONTRASTS"
    label 'process_discovery'   // lightest available label; trivially fast

    // Publish the sentinel only when the threshold is not met so that the
    // orchestrating bash scripts can detect a skipped run without parsing logs.
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true,
               pattern: "low_contrasts.skip"

    input:
    path traitfile
    path boot_traitfile

    output:
    path "traitfile_ok.tab",      emit: traitfile_out,      optional: true
    path "boot_traitfile_ok.tab", emit: boot_traitfile_out, optional: true
    path "low_contrasts.skip",    emit: skip_flag,          optional: true

    script:
    def min_n   = params.min_contrasts ?: 3
    def tname   = params.traitname     ?: 'unknown'
    """
    n_fg=\$(awk '\$2 == 1 { count++ } END { print count+0 }' ${traitfile})

    if [ "\$n_fg" -lt ${min_n} ]; then
        printf "trait=${tname}\\tn_contrasts=\${n_fg}\\tmin_required=${min_n}\\n" \\
            > low_contrasts.skip
        echo "WARNING [CHECK_MIN_CONTRASTS]: Only \${n_fg} foreground contrast(s)" \\
             "in traitfile for trait '${tname}'." \\
             "Minimum ${min_n} required — CT pipeline will be skipped." >&2
    else
        cp ${traitfile}      traitfile_ok.tab
        cp ${boot_traitfile} boot_traitfile_ok.tab
        echo "OK [CHECK_MIN_CONTRASTS]: \${n_fg} foreground contrasts found" \\
             "for trait '${tname}' — proceeding with CT pipeline."
    fi
    """
}
