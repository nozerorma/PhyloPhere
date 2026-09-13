#!/usr/bin/env nextflow
// caas_permulation.nf — Permulation replay for CT discovery significance estimation.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  PERM_REPLAY: Reruns already-permulated labelings (from RESAMPLE) through
 *  caastools' pattern matcher to produce empirical p-value distributions for
 *  CAAS counts.
 *
 *  Consumes:  DISCOVERY output channel (one item per trait/group)
 *  Produces:  perm-replay TSV per trait/group; merged summary fed to CT_CONCAT
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Single-gene perm-replay ────────────────────────────────────────────────────

process PERM_REPLAY {
    label "process_perm_replay"
    tag  "${traitname}/${group}"

    input:
    tuple val(traitname), val(group), path(discovery), path(traitfile), path(tree)

    output:
    tuple val(traitname), val(group), path("perm_replay_${traitname}_${group}.tsv"), emit: perm_replay

    script:
    """
    ct perm-replay \\
        --discovery  ${discovery} \\
        --traitfile  ${traitfile} \\
        --tree       ${tree} \\
        --nboot      ${params.ct_nboot} \\
        --output     perm_replay_${traitname}_${group}.tsv
    """
}


// ── Batched perm-replay (multiple genes per job for efficiency) ───────────────

process PERM_REPLAY_BATCHED {
    label "process_perm_replay_batched"
    tag  "${traitname}/${group}/batch_${batch_id}"

    input:
    tuple val(traitname), val(group), val(batch_id),
          path(discovery_list), path(traitfile), path(tree)

    output:
    tuple val(traitname), val(group), path("perm_replay_batch_${batch_id}_*.tsv"), emit: perm_replay

    script:
    """
    while IFS= read -r disc_file; do
        gene=\$(basename "\${disc_file}" .tsv)
        ct perm-replay \\
            --discovery "\${disc_file}" \\
            --traitfile ${traitfile} \\
            --tree      ${tree} \\
            --nboot     ${params.ct_nboot} \\
            --output    perm_replay_batch_${batch_id}_\${gene}.tsv
    done < ${discovery_list}
    """
}
