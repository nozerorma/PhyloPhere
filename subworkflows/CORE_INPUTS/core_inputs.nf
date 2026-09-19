#!/usr/bin/env nextflow

/*
 * CORE_INPUTS subworkflow
 *
 * Auto-generates the two core reference files the pipeline would otherwise
 * require the user to supply externally, when they are left blank:
 *
 *   - tax_id            (tax_id, species) — resolved from the species tree's
 *                        tip labels via NCBI taxonomy (exact match only).
 *   - gene_ensembl_file (gene, chr, start, end, strand, length,
 *                        human_protein_id) — resolved from the alignment
 *                        gene list via an Ensembl BioMart query.
 *
 * Species/genes that don't resolve exactly are reported (not silently
 * guessed via fuzzy/synonym matching) so the user can supply corrections
 * deliberately. Runs are gated per-file: a param already set by the user is
 * passed through untouched.
 */

process GENERATE_TAXID_MAP {
    tag "auto-generate tax_id map"
    label 'process_low'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path tree_file

    output:
    path "tax_id_generated.tsv",   emit: tax_id_file
    path "tax_id_unresolved.tsv",  emit: unresolved_report

    script:
    """
    python3 ${baseDir}/bin/generate_taxid_map.py \\
        --tree        "${tree_file}" \\
        --output      tax_id_generated.tsv \\
        --unresolved  tax_id_unresolved.tsv
    """
}

process GENERATE_ENSEMBL_MAPPING {
    tag "auto-generate gene_ensembl_file"
    label 'process_low'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path gene_list_file

    output:
    path "gene_ensembl_generated.tsv", emit: ensembl_file
    path "gene_ensembl_unresolved.txt", emit: unresolved_report

    script:
    """
    python3 ${baseDir}/bin/generate_ensembl_mapping.py \\
        --gene-list   "${gene_list_file}" \\
        --output      gene_ensembl_generated.tsv \\
        --unresolved  gene_ensembl_unresolved.txt
    """
}

process DERIVE_GENE_LIST {
    tag "derive gene list from alignment"
    label 'process_low'

    input:
    path alignment_dir

    output:
    path "gene_list.txt", emit: gene_list

    script:
    """
    for f in "${alignment_dir}"/*; do
        [ -f "\$f" ] || continue
        base="\$(basename "\$f")"
        echo "\${base%.*}"
    done | sort -u > gene_list.txt
    """
}

workflow CORE_INPUTS {
    take:
        tree_file   // path/value channel: species tree (newick), or null if unavailable
        alignment_dir_ch // path/value channel: alignment directory, or null if unavailable

    main:
        def tax_id_ch
        def ensembl_ch

        if (params.tax_id) {
            tax_id_ch = Channel.value(file(params.tax_id))
        } else if (tree_file) {
            log.info "[CORE_INPUTS] params.tax_id not set — auto-generating from tree tip labels via NCBI taxonomy."
            GENERATE_TAXID_MAP(tree_file)
            tax_id_ch = GENERATE_TAXID_MAP.out.tax_id_file
        } else {
            log.warn "[CORE_INPUTS] params.tax_id not set and no tree available to auto-generate it from."
            tax_id_ch = Channel.value(file('NO_FILE'))
        }

        if (params.gene_ensembl_file) {
            ensembl_ch = Channel.value(file(params.gene_ensembl_file))
        } else if (alignment_dir_ch) {
            log.info "[CORE_INPUTS] params.gene_ensembl_file not set — auto-generating via Ensembl BioMart."
            DERIVE_GENE_LIST(alignment_dir_ch)
            GENERATE_ENSEMBL_MAPPING(DERIVE_GENE_LIST.out.gene_list)
            ensembl_ch = GENERATE_ENSEMBL_MAPPING.out.ensembl_file
        } else {
            log.warn "[CORE_INPUTS] params.gene_ensembl_file not set and no alignment available to auto-generate it from."
            ensembl_ch = Channel.value(file('NO_FILE'))
        }

    emit:
        tax_id_file          = tax_id_ch
        gene_ensembl_file    = ensembl_ch
}
