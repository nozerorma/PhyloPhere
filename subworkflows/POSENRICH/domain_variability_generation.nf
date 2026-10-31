#!/usr/bin/env nextflow
// domain_variability_generation.nf — Generate the Pfam domain table of POSENRICH when none is supplied.
// PhyloPhere | subworkflows/POSENRICH/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DOMAIN_VARIABILITY_GENERATION: builds the --domain_variability_file by scanning the
 *  reference-species sequence of every gene against Pfam-A with hmmscan
 *  (bin/compute_domain_variability.py). Called from workflows/posenrich.nf.
 *
 *  The Pfam-A.hmm and Pfam-A.clans.tsv cache (about 1.5 GB uncompressed) is downloaded once
 *  to --pfam_cache_dir and reused by later runs.
 *
 *  Consumes:  alignment directory
 *  Produces:  core_inputs/domain_variability.tsv (gene, pfam_id, target_name, description,
 *             clan_acc, clan_name, ali_start, ali_end)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { scoringHash; scoringQuote } from '../SCORING/scoring_contract.nf'


// ── Pfam domain scan ─────────────────────────────────────────────────────────

process COMPUTE_DOMAIN_VARIABILITY {
    cache 'deep'
    tag "auto-generate domain_variability_file"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true,
               pattern: 'domain_variability.tsv*'

    input:
    path alignment_dir

    output:
    path "domain_variability.tsv", emit: domain_variability_file
    path "domain_variability.tsv.manifest.json", emit: manifest, optional: true

    script:
    def cache_dir = params.pfam_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/pfam"
    def source_hash = scoringHash(['compute_domain_variability.py','map_domain_variability.py'].collect { file("${baseDir}/bin/${it}") })
    def prepared = ['Pfam-A.hmm','Pfam-A.clans.tsv','Pfam-A.hmm.h3f','Pfam-A.hmm.h3i','Pfam-A.hmm.h3m','Pfam-A.hmm.h3p'].collect { file("${cache_dir}/${it}") }.findAll { it.isFile() }
    def resource_hash = prepared ? scoringHash(prepared) : 'prepare-on-first-use'
    def ref_species = params.domain_ref_species ?: 'Homo_sapiens'
    """
    # Sources: ${source_hash}; prepared Pfam resource content: ${resource_hash}
    python3 ${baseDir}/bin/compute_domain_variability.py \\
        --alignment-dir "${alignment_dir}" \\
        --output-dir . \\
        --cpus ${task.cpus} \\
        --cache-dir "${cache_dir}" \\
        --ref-species "${ref_species}"
    """

    stub:
    // The script block downloads the Pfam-A cache (about 1.5 GB) on first use, which a stub
    // run must not trigger; without a stub block Nextflow runs the real script even under
    // -stub-run. This stub writes the header of the table only.
    """
    printf 'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\talignment_positions\tcoordinate_space\n' > domain_variability.tsv
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
