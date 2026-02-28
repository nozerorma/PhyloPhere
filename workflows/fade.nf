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
# Tests foreground branches (defined by the caastools traitfile) for
# directional amino-acid selection using HyPhy FADE.
#
# Two directions are always run in parallel:
#   top    → contrast_group == 1 species as foreground (high-trait extreme)
#   bottom → contrast_group == 0 species as foreground (low-trait extreme)
#
# Run modes (params.fade_mode):
#   gene_set  → only genes from CT accumulation / postproc significant lists
#   all       → every PHYLIP alignment file in the alignment directory
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/fade.nf
*/

include { COLLECT_GENE_SETS } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { PHYLIP_TO_FASTA   } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { ANNOTATE_TREE_FG  } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { FADE_RUN              } from "${baseDir}/subworkflows/FADE/fade_run.nf"
include { FADE_REPORT as FADE_REPORT_TOP    } from "${baseDir}/subworkflows/FADE/fade_report.nf"
include { FADE_REPORT as FADE_REPORT_BOTTOM } from "${baseDir}/subworkflows/FADE/fade_report.nf"


// ─── Helpers ─────────────────────────────────────────────────────────────────

/**
 * Build a List of [gene_id, direction, phylip_File] tuples synchronously.
 * Returns a plain Groovy List (not a Channel) so it is safe to call inside
 * a channel operator's closure (e.g. flatMap).
 *
 *   wanted  — Set<String> of gene IDs to keep, or null / empty to keep all.
 */
def ali_tuples_from_dir(String ali_dir, String direction, Set wanted) {
    def ali_path = file(ali_dir)
    if (!ali_path.isDirectory()) {
        log.warn "FADE: alignment directory not found: ${ali_dir} -- skipping '${direction}'"
        return []
    }
    def all_files = ali_path.listFiles()?.findAll { f ->
        f.isFile() && (f.name.endsWith('.phy') ||
                       f.name.endsWith('.phylip') ||
                       f.name.endsWith('.aln') ||
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
        traitfile_input   // Channel<path> or null → falls back to params.caas_config
        tree_input        // Channel<path> or null → falls back to params.tree
        // Optional upstream channels for inline gene-set piping.
        // When non-empty these take priority over the params.fade_* path params.
        // Pass Channel.empty() when running standalone (params-based).
        acc_top_ch        // accumulation CSV for TOP  (*_top_aggregated_results.csv)
        acc_bottom_ch     // accumulation CSV for BOTTOM
        pp_top_ch         // postproc TXT for TOP  (*_change_side_top_significant.txt)
        pp_bottom_ch      // postproc TXT for BOTTOM

    main:

        // ── Resolve shared inputs ────────────────────────────────────────────
        def traitfile_ch = (traitfile_input ?: Channel.empty())
            .ifEmpty {
                assert params.caas_config : "FADE requires a CONTRAST_SELECTION traitfile or --caas_config"
                def f = file(params.caas_config)
                assert f.exists() : "FADE: traitfile not found: ${params.caas_config}"
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

        // ── Build per-direction alignment channels ───────────────────────────
        def all_ali_ch

        if (params.fade_mode == 'all') {

            def top_tuples    = ali_tuples_from_dir(ali_dir, 'top',    null)
            def bottom_tuples = ali_tuples_from_dir(ali_dir, 'bottom', null)
            all_ali_ch = Channel.fromList(top_tuples + bottom_tuples)

        } else {
            // gene_set mode ───────────────────────────────────────────────────
            // Prefer upstream piped channels; fall back to --fade_* path params.
            def resolved_acc_top = (acc_top_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_accumulation_top    ?: 'NO_FILE') }

            def resolved_acc_bottom = (acc_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_accumulation_bottom ?: 'NO_FILE') }

            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_top        ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_bottom     ?: 'NO_FILE') }

            def gene_sets = COLLECT_GENE_SETS(
                resolved_acc_top,
                resolved_acc_bottom,
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

        // ── PHYLIP → FASTA conversion ────────────────────────────────────────
        fasta_ch = PHYLIP_TO_FASTA(all_ali_ch).fasta

        // ── Annotate tree with {Foreground} labels ───────────────────────────
        // fasta_ch is included so that annotate_tree_fg.py can prune the tree
        // to only the taxa present in the per-gene alignment before running FADE
        // (prevents tip-count / sequence-count mismatches in HyPhy FADE).
        def annotate_input_ch = fasta_ch
            .combine(traitfile_ch)
            .combine(tree_ch)
            // → (gene_id, direction, fasta, traitfile, tree)

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

        fade_results_ch = FADE_RUN(fade_input_ch, lg_dat_ch).fade_json

        // ── Reports per direction ────────────────────────────────────────────
        def top_jsons    = fade_results_ch
            .filter { gid, dir, json -> dir == 'top'    }.map { it[2] }.collect().ifEmpty([])
        def bottom_jsons = fade_results_ch
            .filter { gid, dir, json -> dir == 'bottom' }.map { it[2] }.collect().ifEmpty([])

        fade_report_top    = FADE_REPORT_TOP(Channel.value('top'),    top_jsons   )
        fade_report_bottom = FADE_REPORT_BOTTOM(Channel.value('bottom'), bottom_jsons)

    emit:
        report_top     = fade_report_top.report
        report_bottom  = fade_report_bottom.report
        summary_top    = fade_report_top.summary_tsv
        summary_bottom = fade_report_bottom.summary_tsv
        json_results   = fade_results_ch
}
