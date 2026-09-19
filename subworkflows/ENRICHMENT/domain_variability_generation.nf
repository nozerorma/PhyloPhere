#!/usr/bin/env nextflow

/*
 * DOMAIN_VARIABILITY_GENERATION subworkflow
 *
 * Auto-generates ENRICHMENT's --domain_variability_file when left blank, via
 * a cached Pfam-A + hmmscan (bin/compute_domain_variability.py). Reuses
 * ortholog_characterizator's map_domain_variability.py domtblout parser and
 * reference-sequence extraction verbatim (imported, not copied a second
 * time) — the Pfam-A.clans.tsv metadata join producing phylophere's actual
 * schema (gene, pfam_id, target_name, description, clan_acc, clan_name,
 * ali_start, ali_end) is new, since neither repo already builds that file.
 *
 * The Pfam-A.hmm + Pfam-A.clans.tsv cache (~1.5GB uncompressed) is
 * downloaded once to --pfam_cache_dir and reused across runs — never
 * committed ("cache large" pattern, same as STRING/eggNOG).
 */

process COMPUTE_DOMAIN_VARIABILITY {
    tag "auto-generate domain_variability_file"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true,
               pattern: 'domain_variability.tsv'

    input:
    path alignment_dir

    output:
    path "domain_variability.tsv", emit: domain_variability_file

    stub:
    // Explicit stub: the script block downloads a ~1.5GB Pfam-A cache on first
    // use, which -stub-run must never trigger implicitly (no stub: means
    // Nextflow falls back to running the real script block even under
    // -stub-run).
    """
    printf 'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n' > domain_variability.tsv
    """

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
}

workflow DOMAIN_VARIABILITY_GENERATION {
    take:
        alignment_dir_ch  // path/value channel: alignment directory

    main:
        COMPUTE_DOMAIN_VARIABILITY(alignment_dir_ch)

    emit:
        domain_variability_file = COMPUTE_DOMAIN_VARIABILITY.out.domain_variability_file
}
