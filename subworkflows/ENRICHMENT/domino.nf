#!/usr/bin/env nextflow
// domino.nf — Active module identification with DOMINO on a STRING network restricted to the background.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DOMINO_MODULES: finds the active modules of gene lists with DOMINO (Shamir-Lab/DOMINO).
 *  They replace the walktrap clustering of STRING (get_clusters()) and its PPI-density
 *  significance tests (per cluster in describe_clusters(), and for the whole gene list in
 *  run_single_string()), because the significance test of DOMINO accepts the full analysis
 *  background, which get_ppi_enrichment() of STRING does not take at that size. The functional and term enrichment of STRING
 *  (get_enrichment_local() / get_enrichment()) is not replaced and lives in
 *  13.AMI_analysis.Rmd. Called from workflows/enrichment.nf, once per consumer
 *  (caas, fade, rer).
 *
 *  DOMINO_BUILD_NETWORK filters the STRING links file (cached once, unfiltered) to
 *  the edges above domino_network_score_thr whose two endpoints are in the background.
 *  It runs per consumer and not once, because the background differs between consumers
 *  (the universe of RER is not that of FADE).
 *
 *  DOMINO_RUN_MODULES finds and scores the modules of every gene list through the Python
 *  API of DOMINO (run_domino_modules.py calls src.core.domino directly, not the `domino`
 *  CLI), which keeps the Bonferroni-corrected hypergeometric p-value of each module.
 *
 *  Consumes:  background gene list, gene lists, STRING files (string_db_dir or the cache)
 *  Produces:  network.sif, network_edge_scores.tsv, domino_modules/ (published under
 *             ami/domino/<consumer> only with params.publish_domino_intermediates)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Network ──────────────────────────────────────────────────────────────────

process DOMINO_BUILD_NETWORK {
    label 'process_medium'
    publishDir path: "${params.outdir}/ami/domino/${consumer_label}", mode: 'copy', overwrite: true,
               enabled: { params.publish_domino_intermediates ?: false }

    input:
    path background_file
    path gene_list_files  // Not read by the script; it makes the network build wait for the upstream gene lists
    val  score_threshold
    val  string_db_dir
    val  consumer_label   // 'caas' | 'fade' | 'rer'; only names the optional publishDir above

    output:
    path "network.sif",              emit: network_sif
    path "network_edge_scores.tsv",  emit: edge_scores
    path "slices.txt",               emit: slices

    script:
    def db_dir_arg = string_db_dir ? "--string-db-dir ${string_db_dir}" : ""
    // The default --cache-dir of build_domino_network.py ("string_cache") is relative to the
    // work directory of the task, which is new for every task, so the STRING files would be
    // downloaded again on every run. A persistent directory is passed instead.
    def string_cache_dir = params.string_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/string"
    """
    bg_name=\$(basename ${background_file})
    if [ ! -f "${background_file}" ] || [[ "\${bg_name}" == NO_* ]]; then
        echo "[DOMINO_BUILD_NETWORK] Skipping network build: background file is absent or sentinel (\${bg_name})"
        touch network.sif network_edge_scores.tsv slices.txt
    else
        python3 ${baseDir}/subworkflows/ENRICHMENT/local/src/build_domino_network.py \
            --cleaned-background ${background_file} \
            --score-threshold ${score_threshold} \
            ${db_dir_arg} \
            --cache-dir "${string_cache_dir}" \
            --output-dir .

        slicer -n network.sif -o slices.txt
    fi
    """
}

// ── Modules ──────────────────────────────────────────────────────────────────

process DOMINO_RUN_MODULES {
    label 'process_medium'
    publishDir path: "${params.outdir}/ami/domino/${consumer_label}", mode: 'copy', overwrite: true,
               enabled: { params.publish_domino_intermediates ?: false }

    input:
    path network_sif
    path slices_file
    path gene_lists
    val  slice_threshold
    val  module_threshold
    val  consumer_label   // 'caas' | 'fade' | 'rer'; only names the optional publishDir above

    output:
    path "domino_modules", emit: modules_dir

    script:
    """
    mkdir -p domino_gene_lists domino_modules

    if [ ! -s "${network_sif}" ] || [ ! -s "${slices_file}" ]; then
        echo "[DOMINO_RUN_MODULES] Skipping module identification: network or slices file is empty/absent"
        exit 0
    fi

    # The CAAS call site passes a directory of slice_*.tsv files (a header and the gene in column 1)
    if [ -d "${gene_lists}" ]; then
        for f in "${gene_lists}"/slice_*.tsv; do
            [ -f "\$f" ] || continue
            name=\$(basename "\$f" .tsv); name=\${name#slice_}
            tail -n +2 "\$f" | cut -f1 | { grep -v '^[[:space:]]*\$' || true; } > "domino_gene_lists/\${name}.txt"
        done
    fi
    # The FADE and RER call sites pass plain .txt files of active genes, one gene per line
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

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow DOMINO_MODULES {
    take:
    background_file
    gene_list_files
    score_threshold
    slice_threshold
    module_threshold
    string_db_dir
    consumer_label   // 'caas' | 'fade' | 'rer'

    main:
    DOMINO_BUILD_NETWORK(background_file, gene_list_files, score_threshold, string_db_dir, consumer_label)
    DOMINO_RUN_MODULES(
        DOMINO_BUILD_NETWORK.out.network_sif,
        DOMINO_BUILD_NETWORK.out.slices,
        gene_list_files,
        slice_threshold,
        module_threshold,
        consumer_label
    )

    emit:
    network_sif = DOMINO_BUILD_NETWORK.out.network_sif
    edge_scores = DOMINO_BUILD_NETWORK.out.edge_scores
    modules_dir = DOMINO_RUN_MODULES.out.modules_dir
}
