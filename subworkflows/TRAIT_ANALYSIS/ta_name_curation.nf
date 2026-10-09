#!/usr/bin/env nextflow
// ta_name_curation.nf — Rename or prune species-tree tips to match the alignment species names.
// PhyloPhere | subworkflows/TRAIT_ANALYSIS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  NAME_CURATION: curates the input species tree and the trait table so that both use
 *  the species names of the alignment FASTA headers, which they are joined on
 *  downstream. Tips and trait species are translated through a shared NCBI tax_id
 *  (taxonomic synonyms, genus renames); tips with no counterpart in the alignments and
 *  trait species with no counterpart in the tree are removed. The curation is the one
 *  place where synonyms and shared tax_ids are settled: the alphabetically first species
 *  of a shared tax_id keeps it and each other one receives a synthetic tax_id. The
 *  callers (CONTRAST_SELECTION, REPORTING) use the curated tree and trait table in
 *  place of params.tree and params.my_traits.
 *
 *  Species names of the alignments come from:
 *    - params.ali_sp_names, a precomputed flat file (fast path), or
 *    - DERIVE_ALI_SP_NAMES, which scans every FASTA header of params.alignment
 *      (slow path: it reads every alignment file, minutes for thousands of genes).
 *  To avoid the slow path, generate the file once and set params.ali_sp_names:
 *
 *    grep -rh '^>' <alignment_dir> | sed 's/^>//' | sort -u > ali_sp_names.txt
 *
 *  Without alignment names (neither params.ali_sp_names nor params.alignment) every tree
 *  tip is a canonical species and only the trait table is curated.
 *
 *  Consumes:  species tree (Newick), tax_id file (TSV/CSV with tax_id and species,
 *             or a NO_FILE sentinel: names are then matched exactly), trait table
 *             (CSV/TSV, species column params.sp_colname, optional tax_id column) or the
 *             NO_FILE sentinel (only the tree and the species tables are then curated)
 *  Produces:  curated_tree (Newick); report (TSV: original_name, curated_name,
 *             fate = kept, renamed or pruned, per tip); curated_traits (the trait table
 *             with species renamed to the tip names and unmatched rows removed);
 *             species_table (TSV: species, tax_id, tax_id_resolved, note);
 *             taxid_map (TSV in the layout of params.tax_id: one row per canonical species
 *             with its resolved tax_id and family);
 *             species_report (text: source, original_name, curated_name, status =
 *             maintained, changed or removed, and the reason), published in name_curation/
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
    path trait_file

    output:
    path "curated_tree.nwk",                 emit: curated_tree
    path "name_curation_report.tsv",         emit: report
    path "curated_traits.*",                 emit: curated_traits, optional: true
    path "species_table.tsv",                emit: species_table
    path "species_curation_report.txt",      emit: species_report
    path "species_taxid_map.tsv",            emit: taxid_map

    script:
    // The curated trait table keeps the format (extension) of the input table. Without a trait
    // table (the NO_FILE sentinel) only the tree and the species tables are curated.
    def with_traits = trait_file.name != 'NO_FILE'
    def ext = trait_file.name.tokenize('.').last()
    def trait_args = with_traits
        ? "--traits \"${trait_file}\" --sp-col \"${params.sp_colname}\" --traits-out \"curated_traits.${ext}\""
        : ""
    """
    python3 ${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/tree_cleanup.py \\
        --tree          "${tree_file}" \\
        --ali-sp-names  "${ali_sp_names_file}" \\
        --tax-id        "${tax_id_file}" \\
        --output        curated_tree.nwk \\
        --report        name_curation_report.tsv \\
        ${trait_args} \\
        --species-table species_table.tsv \\
        --species-report species_curation_report.txt \\
        --taxid-map     species_taxid_map.tsv
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow NAME_CURATION {
    take:
        tree_file    // path channel: input species tree (newick)
        tax_id_file  // path/value channel: taxid-to-species TSV
        trait_file   // path/value channel: trait table (CSV/TSV), or the NO_FILE sentinel

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
            log.warn "[NAME_CURATION] Neither params.ali_sp_names nor params.alignment is set: " +
                     "the tree tips are taken as canonical species and only the trait table is curated."
            ali_sp_ch = Channel.value(file('NO_FILE'))
        }

        TREE_CLEANUP(tree_file, tax_id_file, ali_sp_ch, trait_file)

    emit:
        curated_tree   = TREE_CLEANUP.out.curated_tree
        report         = TREE_CLEANUP.out.report
        curated_traits = TREE_CLEANUP.out.curated_traits
        species_table  = TREE_CLEANUP.out.species_table
        species_report = TREE_CLEANUP.out.species_report
        taxid_map      = TREE_CLEANUP.out.taxid_map
}
