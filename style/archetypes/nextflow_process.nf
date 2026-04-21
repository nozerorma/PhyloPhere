#!/usr/bin/env nextflow
// ct_bootstrap.nf — Bootstrap resampling for CT discovery significance estimation.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_BOOTSTRAP: Runs caastools bootstrap on per-group discovery outputs to produce
 *  empirical p-value distributions for CAAS counts.
 *
 *  Consumes:  DISCOVERY output channel (one item per trait/group)
 *  Produces:  bootstrap TSV per trait/group; merged summary fed to CT_CONCAT
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Single-gene bootstrap ─────────────────────────────────────────────────────

process BOOTSTRAP {
    label "process_boot"
    tag  "${traitname}/${group}"

    input:
    tuple val(traitname), val(group), path(discovery), path(traitfile), path(tree)

    output:
    tuple val(traitname), val(group), path("boot_${traitname}_${group}.tsv"), emit: bootstrap

    script:
    """
    caastools bootstrap \\
        --discovery  ${discovery} \\
        --traitfile  ${traitfile} \\
        --tree       ${tree} \\
        --nboot      ${params.ct_nboot} \\
        --output     boot_${traitname}_${group}.tsv
    """
}


// ── Batched bootstrap (multiple genes per job for efficiency) ─────────────────

process BOOTSTRAP_BATCHED {
    label "process_boot_batched"
    tag  "${traitname}/${group}/batch_${batch_id}"

    input:
    tuple val(traitname), val(group), val(batch_id),
          path(discovery_list), path(traitfile), path(tree)

    output:
    tuple val(traitname), val(group), path("boot_batch_${batch_id}_*.tsv"), emit: bootstrap

    script:
    """
    while IFS= read -r disc_file; do
        gene=\$(basename "\${disc_file}" .tsv)
        caastools bootstrap \\
            --discovery "\${disc_file}" \\
            --traitfile ${traitfile} \\
            --tree      ${tree} \\
            --nboot     ${params.ct_nboot} \\
            --output    boot_batch_${batch_id}_\${gene}.tsv
    done < ${discovery_list}
    """
}
