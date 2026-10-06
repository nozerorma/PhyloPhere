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

    // The sentinel exists only when the threshold is not met, so a launch script can detect a
    // skipped run without parsing logs.
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true,
               pattern: "low_contrasts.skip"

    input:
    path traitfile
    path permulation_traitfile
    path trait_dir, stageAs: 'trait_dir_in'

    output:
    path "traitfile_ok.tab",             emit: traitfile_out,             optional: true
    path "permulation_traitfile_ok.tab", emit: permulation_traitfile_out, optional: true
    path "traitfiles_ok_dir",            emit: trait_dir_out,             optional: true
    path "low_contrasts.skip",           emit: skip_flag,                 optional: true

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
        cp ${traitfile}             traitfile_ok.tab
        cp ${permulation_traitfile} permulation_traitfile_ok.tab
        mkdir -p traitfiles_ok_dir
        # params.multi_hypothesis chooses between all the traitfile_H*.tab hypotheses and the canonical one (H1)
        if [ "${params.multi_hypothesis}" = "true" ] && [ -d "${trait_dir}" ]; then
            cp ${trait_dir}/traitfile_H*.tab traitfiles_ok_dir/ 2>/dev/null || cp ${traitfile} traitfiles_ok_dir/traitfile_H1.tab
        else
            if [ -f "${trait_dir}/traitfile_H1.tab" ]; then
                cp ${trait_dir}/traitfile_H1.tab traitfiles_ok_dir/traitfile_H1.tab
            else
                cp ${traitfile} traitfiles_ok_dir/traitfile_H1.tab
            fi
        fi
        # The per-pair PSS weights (contrast_hypotheses_pairs.tsv, written by
        # 4.Independent_contrasts.Rmd) travel with the passing traitfiles: traitfiles_ok_dir
        # is the source of the domain-pool weights of SCORING (scoring_hyp_pairs_ch in main.nf).
        # The copy is allowed to fail when the files do not exist.
        if [ -d "${trait_dir}" ]; then
            cp ${trait_dir}/contrast_hypotheses_pairs.tsv traitfiles_ok_dir/ 2>/dev/null || true
            cp ${trait_dir}/contrast_hypotheses_summary.tsv traitfiles_ok_dir/ 2>/dev/null || true
        fi
        echo "OK [CHECK_MIN_CONTRASTS]: \${n_fg} foreground contrasts found" \\
             "for trait '${tname}' — proceeding with CT pipeline."
    fi
    """

}
