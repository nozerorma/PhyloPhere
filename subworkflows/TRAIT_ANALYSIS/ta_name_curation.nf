#!/usr/bin/env nextflow
// ta_name_curation.nf — Rename or prune species-tree tips to match the alignment species names.
// PhyloPhere | subworkflows/TRAIT_ANALYSIS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  NAME_CURATION: curates the input species tree so that its tip labels match the
 *  species names of the alignment FASTA headers, which the traitfile and the
 *  alignments are joined on downstream. Tips are translated through a shared NCBI
 *  tax_id (taxonomic synonyms, genus renames) and tips with no counterpart in the
 *  alignments are pruned. The callers (CONTRAST_SELECTION, REPORTING) use the
 *  curated tree in place of params.tree.
 *
 *  Species names of the alignments come from:
 *    - params.ali_sp_names, a precomputed flat file (fast path), or
 *    - DERIVE_ALI_SP_NAMES, which scans every FASTA header of params.alignment
 *      (slow path: it reads every alignment file, minutes for thousands of genes).
 *  To avoid the slow path, generate the file once and set params.ali_sp_names:
 *
 *    grep -rh '^>' <alignment_dir> | sed 's/^>//' | sort -u > ali_sp_names.txt
 *
 *  Consumes:  species tree (Newick), tax_id file (TSV/CSV with tax_id and species,
 *             or a NO_FILE sentinel: tips are then matched by exact name only)
 *  Produces:  curated_tree (Newick) and report (TSV: original_name, curated_name,
 *             fate = kept, renamed or pruned, per tip), published in name_curation/
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// ── Processes ────────────────────────────────────────────────────────────────

process DERIVE_ALI_SP_NAMES {
    tag "derive species names from alignments"
    label 'process_medium'

    input:
    path alignment_dir

    output:
    path "ali_sp_names.txt", emit: sp_names

    script:
    """
    echo "[NAME_CURATION] Deriving species names by scanning alignment directory."
    echo "[NAME_CURATION] WARNING: this scans all alignment files and may be slow."
    grep -rh "^>" "${alignment_dir}" 2>/dev/null \
        | sed 's/^>//' \
        | tr -d '\\r' \
        | sort -u \
        > ali_sp_names.txt
    n=\$(wc -l < ali_sp_names.txt | tr -d ' ')
    echo "[NAME_CURATION] Found \${n} unique species names in alignment directory."
    """
}

process TREE_CLEANUP {
    tag "curate tree tip labels"
    label 'process_low'

    publishDir "${params.outdir}/name_curation", mode: 'copy', overwrite: true

    input:
    path tree_file
    path tax_id_file
    path ali_sp_names_file

    output:
    path "curated_tree.nwk",         emit: curated_tree
    path "name_curation_report.tsv", emit: report

    script:
    """
    python3 ${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/tree_cleanup.py \\
        --tree          "${tree_file}" \\
        --ali-sp-names  "${ali_sp_names_file}" \\
        --tax-id        "${tax_id_file}" \\
        --output        curated_tree.nwk \\
        --report        name_curation_report.tsv
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow NAME_CURATION {
    take:
        tree_file    // path channel: input species tree (newick)
        tax_id_file  // path/value channel: taxid-to-species TSV

    main:
        def ali_sp_ch

        if (params.ali_sp_names) {
            log.info "[NAME_CURATION] Using pre-built species list: ${params.ali_sp_names}"
            ali_sp_ch = Channel.value(file(params.ali_sp_names))
        } else if (params.alignment) {
            log.warn """\
                [NAME_CURATION] params.ali_sp_names not set — deriving species names by scanning
                the alignment directory (${params.alignment}).
                This reads every alignment file and may be very slow for large datasets.
                To avoid this cost on future runs, generate the file once:
                  grep -rh '^>' ${params.alignment} | sed 's/^>//' | sort -u > ali_sp_names.txt
                then set params.ali_sp_names in your config.
                """.stripIndent()
            DERIVE_ALI_SP_NAMES(Channel.value(file(params.alignment, type: 'dir')))
            ali_sp_ch = DERIVE_ALI_SP_NAMES.out.sp_names
        } else {
            error "[NAME_CURATION] Neither params.ali_sp_names nor params.alignment is set. " +
                  "Provide at least one to run name curation."
        }

        TREE_CLEANUP(tree_file, tax_id_file, ali_sp_ch)

    emit:
        curated_tree = TREE_CLEANUP.out.curated_tree
        report       = TREE_CLEANUP.out.report
}
