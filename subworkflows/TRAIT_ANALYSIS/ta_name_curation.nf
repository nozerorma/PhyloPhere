#!/usr/bin/env nextflow

/*
 * NAME_CURATION subworkflow
 *
 * Curates the input species tree so that all tip labels match the naming
 * convention used in the alignment files. This is a prerequisite for
 * downstream CT, RERConverge, and FADE analyses that join tree-derived
 * species names (traitfile) with alignment FASTA headers.
 *
 * Strategy
 * --------
 * 1. Obtain the canonical set of species names present in the alignments:
 *      a. If params.ali_sp_names is provided: use it directly (fast path).
 *      b. Otherwise: scan all FASTA files in the alignment directory to
 *         extract unique header labels (slow path — reads every alignment
 *         file; see warning below).
 * 2. Use params.tax_id as a shared NCBI taxon-ID key to translate tree tip
 *    labels that differ from alignment names (e.g. taxonomic synonyms).
 * 3. Prune tree tips that have no match in the alignment species set after
 *    translation.
 *
 * Output
 * ------
 *   curated_tree  —  Newick tree whose tip labels exactly match alignment
 *                    headers; passed downstream as the canonical tree,
 *                    replacing params.tree everywhere.
 *   report        —  TSV summarising each tip: kept / renamed / pruned.
 *
 * WARNING (slow path)
 * -------------------
 * When ali_sp_names is NOT provided, this subworkflow scans every alignment
 * file in the alignment directory to build the species list. For large
 * datasets (thousands of genes × hundreds of species) this can take several
 * minutes. It is strongly recommended to pre-generate the file once with:
 *
 *   grep -rh '^>' <alignment_dir> | sed 's/^>//' | sort -u > ali_sp_names.txt
 *
 * and set params.ali_sp_names in your config to avoid this overhead on every
 * run.
 */

// ─── Processes ────────────────────────────────────────────────────────────────

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
    python3 ${baseDir}/subworkflows/TRAIT_ANALYSIS/local/scripts/tree_cleanup.py \\
        --tree          "${tree_file}" \\
        --ali-sp-names  "${ali_sp_names_file}" \\
        --tax-id        "${tax_id_file}" \\
        --output        curated_tree.nwk \\
        --report        name_curation_report.tsv
    """
}

// ─── Workflow ─────────────────────────────────────────────────────────────────

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
