#!/usr/bin/env nextflow
// domain_variability_generation.nf — Generate the Pfam domain table of ENRICHMENT when none is supplied.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DOMAIN_VARIABILITY_GENERATION: builds the --domain_variability_file by scanning the
 *  reference-species sequence of every gene against Pfam-A with hmmscan
 *  (bin/compute_domain_variability.py). Called from workflows/enrichment.nf.
 *
 *  The Pfam-A.hmm and Pfam-A.clans.tsv cache (about 1.5 GB uncompressed) is downloaded once
 *  to --pfam_cache_dir and reused by later runs.
 *
 *  Consumes:  alignment directory
 *  Produces:  core_inputs/domain_variability.tsv (gene, pfam_id, target_name, description,
 *             clan_acc, clan_name, ali_start, ali_end)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Pfam domain scan ─────────────────────────────────────────────────────────

process COMPUTE_DOMAIN_VARIABILITY {
    tag "auto-generate domain_variability_file"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true,
               pattern: 'domain_variability.tsv'

    input:
    path alignment_dir

    output:
    path "domain_variability.tsv", emit: domain_variability_file

    script:
    def cache_dir = params.pfam_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/pfam"
    def ref_species = params.domain_ref_species ?: 'Homo_sapiens'
    """
    python3 ${baseDir}/bin/compute_domain_variability.py \\
        --alignment-dir "${alignment_dir}" \\
        --output-dir . \\
        --cache-dir "${cache_dir}" \\
        --ref-species "${ref_species}"
    """

    stub:
    // The script block downloads the Pfam-A cache (about 1.5 GB) on first use, which a stub
    // run must not trigger; without a stub block Nextflow runs the real script even under
    // -stub-run. This stub writes the header of the table only.
    """
    printf 'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n' > domain_variability.tsv
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow DOMAIN_VARIABILITY_GENERATION {
    take:
        alignment_dir_ch  // path/value channel: alignment directory

    main:
        COMPUTE_DOMAIN_VARIABILITY(alignment_dir_ch)

    emit:
        domain_variability_file = COMPUTE_DOMAIN_VARIABILITY.out.domain_variability_file
}
