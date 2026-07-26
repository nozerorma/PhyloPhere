#!/usr/bin/env nextflow

/*
 * DOMINO — active module identification
 * ────────────────────────────────────────────────────────────────────
 * Replaces STRING's walktrap clustering (get_clusters()) and both its
 * PPI-density significance tests (per-cluster inside describe_clusters(), and
 * the whole-gene-list test in run_single_string()) with DOMINO
 * (Shamir-Lab/DOMINO): its module-finding recovers known complexes more fully
 * than walktrap at our list sizes (10-40 genes), and its significance test
 * accepts an arbitrary background where STRING's get_ppi_enrichment() cannot
 * (STRING hard-caps foreground+background at 2000 proteins combined; our
 * cleaned_background is ~10,447). STRING's functional/term enrichment
 * (get_enrichment_local()/get_enrichment()) is untouched and lives entirely in
 * 15.AMI_analysis.Rmd.
 *
 * DOMINO_BUILD_NETWORK filters STRING v12.0's raw links file (cached once,
 * unfiltered) to domino_network_score_thr and both-endpoints-in-background.
 * This has to be rebuilt per report, not shared globally, because the
 * background varies by consumer (e.g. RER's universe differs from FADE's —
 * see bugfix_fcs_universe_wrong_background.md).
 *
 * DOMINO_RUN_MODULES finds + scores modules per gene list via DOMINO's own
 * Python API (run_domino_modules.py calls src.core.domino directly, not the
 * `domino` CLI, so the Bonferroni-corrected hypergeometric p-value per module
 * survives instead of being discarded the way modules.out is).
 */

process DOMINO_BUILD_NETWORK {
    label 'process_medium'

    input:
    path background_file
    val  score_threshold
    val  string_db_dir

    output:
    path "network.sif", emit: network_sif
    path "slices.txt",  emit: slices

    script:
    def db_dir_arg = string_db_dir ? "--string-db-dir ${string_db_dir}" : ""
    """
    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/build_domino_network.py \
        --cleaned-background ${background_file} \
        --score-threshold ${score_threshold} \
        ${db_dir_arg} \
        --output-dir .

    slicer -n network.sif -o slices.txt
    """
}

process DOMINO_RUN_MODULES {
    label 'process_medium'

    input:
    path network_sif
    path slices_file
    path gene_lists
    val  slice_threshold
    val  module_threshold

    output:
    path "domino_modules", emit: modules_dir

    script:
    """
    mkdir -p domino_gene_lists domino_modules

    # SCORING call site: gene_lists is a directory of slice_*.tsv (Gene column + header)
    if [ -d "${gene_lists}" ]; then
        for f in "${gene_lists}"/slice_*.tsv; do
            [ -f "\$f" ] || continue
            name=\$(basename "\$f" .tsv); name=\${name#slice_}
            tail -n +2 "\$f" | cut -f1 | { grep -v '^[[:space:]]*\$' || true; } > "domino_gene_lists/\${name}.txt"
        done
    fi
    # TOOL_* call sites: gene_lists are already flat, plain active-gene .txt files
    for f in *.txt; do
        [ -f "\$f" ] || continue
        cp "\$f" domino_gene_lists/ 2>/dev/null || true
    done

    python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/run_domino_modules.py \
        --network ${network_sif} \
        --slices ${slices_file} \
        --gene-lists-dir domino_gene_lists \
        --slice-threshold ${slice_threshold} \
        --module-threshold ${module_threshold} \
        --threads ${task.cpus} \
        --output-dir domino_modules
    """
}

workflow DOMINO_MODULES {
    take:
    background_file
    gene_list_files
    score_threshold
    slice_threshold
    module_threshold
    string_db_dir

    main:
    DOMINO_BUILD_NETWORK(background_file, score_threshold, string_db_dir)
    DOMINO_RUN_MODULES(
        DOMINO_BUILD_NETWORK.out.network_sif,
        DOMINO_BUILD_NETWORK.out.slices,
        gene_list_files,
        slice_threshold,
        module_threshold
    )

    emit:
    network_sif = DOMINO_BUILD_NETWORK.out.network_sif
    modules_dir = DOMINO_RUN_MODULES.out.modules_dir
}
