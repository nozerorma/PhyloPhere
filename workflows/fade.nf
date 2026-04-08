#!/usr/bin/env nextflow

/*
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: FADE Workflow
#
# Tests foreground branches (defined by trait_stats.csv extremes) for
# directional amino-acid selection using HyPhy FADE.
#
# Two directions are always run in parallel:
#   top    → global_label == high_extreme species as foreground
#   bottom → global_label == low_extreme species as foreground
#
# Run modes (params.fade_mode):
#   gene_set  → only genes from CT postproc significant lists
#   all       → every supported alignment file in the alignment directory
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/fade.nf
*/

include { COLLECT_GENE_SETS       } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { EXTRACT_EXTREME_SPECIES } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { PHYLIP_TO_FASTA         } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { FILTER_FASTA_TO_TREE    } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { ANNOTATE_TREE_FG        } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { FADE_RUN; FADE_BATCHED } from "${baseDir}/subworkflows/FADE/fade_run.nf"
include { FADE_REPORT as FADE_REPORT_TOP    } from "${baseDir}/subworkflows/FADE/fade_report.nf"
include { FADE_REPORT as FADE_REPORT_BOTTOM } from "${baseDir}/subworkflows/FADE/fade_report.nf"


// ─── Helpers ─────────────────────────────────────────────────────────────────

/**
 * Render a list of TSV row strings into a heredoc-safe manifest block.
 * Strips leading whitespace (Groovy multiline strings indent with the closure).
 */
def createBatchManifestText = { List<String> rows ->
    rows
        .collect { row -> row.replaceFirst(/^\s+/, '') }
        .join(System.lineSeparator()) + System.lineSeparator()
}

/**
 * Build a List of [gene_id, direction, alignment_File] tuples synchronously.
 * Returns a plain Groovy List (not a Channel) so it is safe to call inside
 * a channel operator's closure (e.g. flatMap).
 *
 *   wanted  — Set<String> of gene IDs to keep, or null / empty to keep all.
 */
def ali_tuples_from_dir(String ali_dir, String direction, Set wanted) {
    def extensions = ['.phy', '.phylip', '.aln', '.fa', '.fasta']

    // ── plain directory ─────────────────────────────────────────────────────
    def ali_path = file(ali_dir)
    if (!ali_path.isDirectory()) {
        log.warn "FADE: alignment directory not found: ${ali_dir} -- skipping '${direction}'"
        return []
    }
    def all_files = ali_path.listFiles()?.findAll { f ->
        def name = f.name.toLowerCase()
        f.isFile() && (name.endsWith('.phy') ||
                       name.endsWith('.phylip') ||
                       name.endsWith('.aln') ||
                       name.endsWith('.fa') ||
                       name.endsWith('.fasta') ||
                       !f.name.contains('.'))
    } ?: []

    def tuples = all_files.collect { f ->
        def gid = f.name.tokenize('.')[0]
        (wanted == null || wanted.contains(gid)) ? [gid, direction, f] : null
    }.findAll { it != null }

    if (tuples.isEmpty())
        log.warn "FADE: no alignment files matched for direction '${direction}' in ${ali_dir}"
    else
        log.info "FADE: ${tuples.size()} genes queued for '${direction}'"
    return tuples
}


// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow FADE {

    take:
        stats_file_input  // Channel<path> or null → falls back to canonical trait_stats.csv under outdir
        tree_input        // Channel<path> or null → falls back to params.tree
        // Optional upstream channels for inline gene-set piping.
        // When non-empty these take priority over the params.fade_* path params.
        // Pass Channel.empty() when running standalone (params-based).
        pp_top_ch         // postproc TXT for TOP  (all_top.txt)
        pp_bottom_ch      // postproc TXT for BOTTOM (all_bottom.txt)
        // Optional CT discovery table used to reuse the single upstream toy sample.
        ct_discovery_input

    main:

        // ── Resolve shared inputs ────────────────────────────────────────────
        def stats_file_ch = (stats_file_input ?: Channel.empty())
            .ifEmpty {
                def fallback = "${params.outdir}/data_exploration/1.Data-exploration/1.Species_distribution/trait_stats.csv"
                def f = file(fallback)
                assert f.exists() : "FADE: trait_stats.csv not found at ${fallback}"
                f
            }

        def tree_ch = (tree_input ?: Channel.empty())
            .ifEmpty {
                assert params.tree : "FADE requires --tree"
                def f = file(params.tree)
                assert f.exists() : "FADE: tree not found: ${params.tree}"
                f
            }

        def ali_dir = (params.fade_alignment ?: params.alignment) ?:
            error("FADE: no alignment directory specified (--fade_alignment or --alignment)")

        // ── LG model dat file ────────────────────────────────────────────────
        def lg_dat_ch = Channel.value(file(params.lg_dat_path))
        def species_lists = EXTRACT_EXTREME_SPECIES(stats_file_ch)
        def top_species_ch = species_lists.top_species
        def bottom_species_ch = species_lists.bottom_species

        // ── Build per-direction alignment channels ───────────────────────────
        def all_ali_ch

        if (params.fade_mode == 'all') {
            if (params.toy_mode) {
                def resolved_ct_discovery = (ct_discovery_input ?: Channel.empty())
                    .ifEmpty { file('NO_FILE') }

                all_ali_ch = resolved_ct_discovery.flatMap { discovery_file ->
                    def toy_genes = null
                    if (discovery_file.name != 'NO_FILE' && discovery_file.exists()) {
                        toy_genes = discovery_file
                            .readLines()
                            .drop(1)
                            .collect { it.split('\t')[0].trim() }
                            .findAll { it }
                            .toSet()
                        log.info "[toy_mode] FADE: reusing ${toy_genes.size()} genes from CT discovery"
                    } else {
                        def n = (params.toy_n ?: 50) as int
                        def all_top = ali_tuples_from_dir(ali_dir, 'top', null)
                        Collections.shuffle(all_top)
                        toy_genes = all_top.take(n).collect { it[0] }.toSet()
                        log.info "[toy_mode] FADE: using ${toy_genes.size()} randomly sampled genes"
                    }

                    def top_tuples    = ali_tuples_from_dir(ali_dir, 'top',    toy_genes)
                    def bottom_tuples = ali_tuples_from_dir(ali_dir, 'bottom', toy_genes)
                    top_tuples + bottom_tuples
                }
            } else {
                def top_tuples    = ali_tuples_from_dir(ali_dir, 'top',    null)
                def bottom_tuples = ali_tuples_from_dir(ali_dir, 'bottom', null)
                all_ali_ch = Channel.fromList(top_tuples + bottom_tuples)
            }

        } else {
            // gene_set mode ───────────────────────────────────────────────────
            // Prefer upstream piped channels; fall back to --fade_* path params.
            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_top        ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_bottom     ?: 'NO_FILE') }

            def gene_sets = COLLECT_GENE_SETS(
                resolved_pp_top,
                resolved_pp_bottom
            )

            // flatMap runs the closure synchronously at evaluation time, where
            // file I/O is valid. Returns a List, not a Channel, so flatMap
            // correctly emits each tuple as a separate item.
            def top_ali_ch = gene_sets.gene_set_top.flatMap { gsf ->
                def wanted = gsf.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#') }.toSet()
                ali_tuples_from_dir(ali_dir, 'top', wanted)
            }

            def bottom_ali_ch = gene_sets.gene_set_bottom.flatMap { gsf ->
                def wanted = gsf.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#') }.toSet()
                ali_tuples_from_dir(ali_dir, 'bottom', wanted)
            }

            all_ali_ch = top_ali_ch.mix(bottom_ali_ch)
        }

        // ── Normalize all supported alignment inputs to FASTA ────────────────
        // FASTA files (.fa / .fasta) bypass the conversion process entirely;
        // only PHYLIP-like inputs are submitted to PHYLIP_TO_FASTA.
        def ali_branched = all_ali_ch.branch {
            is_fasta: { gid, dir, f ->
                def n = f.name.toLowerCase()
                n.endsWith('.fa') || n.endsWith('.fasta')
            }
            needs_convert: true
        }
        def direct_fasta_ch  = ali_branched.is_fasta
        def converted_ch     = PHYLIP_TO_FASTA(ali_branched.needs_convert).fasta
        fasta_ch = direct_fasta_ch.mix(converted_ch)

        // Restrict every gene alignment to the phenotype-pruned tree taxa
        // before FADE-specific foreground labeling happens.
        def fasta_in_tree_ch = fasta_ch
            .combine(tree_ch)
            .map { gid, dir, fa, tree -> tuple(gid, dir, fa, tree) }
        def tree_filtered_fasta_ch = FILTER_FASTA_TO_TREE(fasta_in_tree_ch).fasta

        // ── Annotate tree with {Foreground} labels ───────────────────────────
        // fasta_ch is included so that annotate_tree_fg.py can prune the tree
        // to only the taxa present in the per-gene alignment before running FADE
        // (prevents tip-count / sequence-count mismatches in HyPhy FADE).
        def top_annotate_input_ch = tree_filtered_fasta_ch
            .filter { gid, dir, fa -> dir == 'top' }
            .combine(top_species_ch)
            .combine(tree_ch)

        def bottom_annotate_input_ch = tree_filtered_fasta_ch
            .filter { gid, dir, fa -> dir == 'bottom' }
            .combine(bottom_species_ch)
            .combine(tree_ch)

        def annotate_input_ch = top_annotate_input_ch.mix(bottom_annotate_input_ch)
            // → (gene_id, direction, fasta, fg_species_file, tree)

        def annotate_result    = ANNOTATE_TREE_FG(annotate_input_ch)
        annotated_ch           = annotate_result.annotated_tree
        def filtered_fasta_ch  = annotate_result.filtered_fasta

        // ── Join filtered FASTA + annotated tree, then run FADE ─────────────
        //  filtered_fasta_ch : (gene_id, direction, filtered_fasta)
        //    — FASTA pruned to only sequences present in the (pruned) tree;
        //      prevents HyPhy tip-count / sequence-count mismatch in both
        //      directions (tree > fasta AND fasta > tree).
        //  annotated_ch      : (gene_id, direction, annotated_tree)
        //  join by [0,1]     → (gene_id, direction, filtered_fasta, annotated_tree)
        def fade_input_ch = filtered_fasta_ch.join(annotated_ch, by: [0, 1])

        def fadeBatchSize = (params.fade_batch_size ?: 1) as int
        if (fadeBatchSize > 1) {
            // ── Batched mode ──────────────────────────────────────────────
            // Split by direction first so every batch is mono-directional,
            // then collate each direction's genes into chunks of fadeBatchSize.
            def fade_branched = fade_input_ch.branch {
                top:    it[1] == 'top'
                bottom: it[1] == 'bottom'
            }
            def fadeBatchCounterTop    = 0
            def fadeBatchCounterBottom = 0

            def make_fade_batches = { branch_ch, dir, counter ->
                branch_ch
                    .collate(fadeBatchSize)
                    .map { batch ->
                        def batchID = sprintf("fade_batch_${dir}_%05d", ++counter)
                        def manifestText = createBatchManifestText(
                            batch.collect { row -> "${row[0]}\t${row[2].name}\t${row[3].name}" }
                        )
                        tuple(batchID, dir, batch.size(), manifestText,
                              batch.collect { row -> row[2] },   // filtered fastas
                              batch.collect { row -> row[3] })   // annotated trees
                    }
            }

            def batches_ch = make_fade_batches(fade_branched.top,    'top',    fadeBatchCounterTop)
                .mix(make_fade_batches(fade_branched.bottom, 'bottom', fadeBatchCounterBottom))

            def batched_out = FADE_BATCHED(batches_ch, lg_dat_ch).fade_json
            // fade_json emits (direction, path) per JSON; flatten to (gene_id, direction, path)
            fade_results_ch = batched_out
                .transpose()
                .map { dir, f ->
                    def gene_id = f.name.replace(".${dir}.FADE.json", "")
                    tuple(gene_id, dir, f)
                }
        } else {
            // ── Single-gene mode (default) ────────────────────────────────
            fade_results_ch = FADE_RUN(fade_input_ch, lg_dat_ch).fade_json
        }

        // ── Reports per direction ────────────────────────────────────────────
        //    Use .branch() to split once and guarantee exclusive routing,
        //    avoiding any ambiguity from applying two independent .filter()
        //    subscriptions on the same queue channel.
        def branched = fade_results_ch.branch {
            top:    it[1] == 'top'
            bottom: it[1] == 'bottom'
        }
        def top_jsons    = branched.top.map    { it[2] }.collect().ifEmpty([])
        def bottom_jsons = branched.bottom.map { it[2] }.collect().ifEmpty([])

        fade_report_top    = FADE_REPORT_TOP(Channel.value('top'),    top_jsons   )
        fade_report_bottom = FADE_REPORT_BOTTOM(Channel.value('bottom'), bottom_jsons)

    emit:
        report_top     = fade_report_top.report
        report_bottom  = fade_report_bottom.report
        summary_top    = fade_report_top.summary_tsv
        summary_bottom = fade_report_bottom.summary_tsv
        json_results   = fade_results_ch
}
