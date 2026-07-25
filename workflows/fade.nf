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
# Alignment preparation (PHYLIP→FASTA conversion, tree-filtering, and species
# extraction) is handled upstream by SELECTION_PREP (called once in main.nf)
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/fade.nf
*/

include { ANNOTATE_TREE_FG; ANNOTATE_TREE_FG_BATCHED } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { FADE_RUN; FADE_BATCHED  } from "${baseDir}/subworkflows/FADE/fade_run.nf"
include { FADE_REPORT as FADE_REPORT_TOP    } from "${baseDir}/subworkflows/FADE/fade_report.nf"
include { FADE_REPORT as FADE_REPORT_BOTTOM } from "${baseDir}/subworkflows/FADE/fade_report.nf"
include { FADE_GENE_LISTS as FADE_GENE_LISTS_TOP; FADE_GENE_LISTS as FADE_GENE_LISTS_BOTTOM } from "${baseDir}/subworkflows/FADE/fade_gene_lists.nf"
include { FADE_JSON_TO_CSV as FADE_JSON_TO_CSV_TOP; FADE_JSON_TO_CSV as FADE_JSON_TO_CSV_BOTTOM } from "${baseDir}/subworkflows/FADE/fade_json_to_csv.nf"
include { TOOL_STRING_MODULE_REPORT as FADE_STRING_REPORT } from "${baseDir}/subworkflows/ENRICHMENT/tool_enrichment.nf"


// ─── Helpers ─────────────────────────────────────────────────────────────────

/**
 * Render a list of TSV row strings into a heredoc-safe manifest block.
 */
def createBatchManifestText = { List<String> rows ->
    rows
        .collect { row -> row.replaceFirst(/^\s+/, '') }
        .join(System.lineSeparator()) + System.lineSeparator()
}


// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow FADE {

    take:
        // Pre-processed alignment channels from SELECTION_PREP (called once in main.nf).
        // Each channel emits (gene_id, filtered_fasta) tuples — no direction label yet;
        // direction is introduced inside this workflow for ANNOTATE_TREE_FG.
        fasta_top_ch      // (gene_id, fasta) for top-direction genes
        fasta_bottom_ch   // (gene_id, fasta) for bottom-direction genes
        // Shared reference files (value channels — can be broadcast to every gene)
        tree_ch           // path to phenotype-pruned species tree
        top_species_ch    // path to top_species.txt
        bottom_species_ch // path to bottom_species.txt

    main:

        // ── LG model dat file ────────────────────────────────────────────────
        def lg_dat_ch = Channel.value(file(params.lg_dat_path))

        // ── Annotate tree with {Foreground} labels ───────────────────────────
        // Introduce direction labels and combine with the species/tree files.
        // ANNOTATE_TREE_FG also prunes the tree and FASTA to mutual taxa,
        // preventing HyPhy tip-count / sequence-count mismatches.
        //
        // Input tuple: (gene_id, direction, fasta, fg_species_file, tree)
        def top_annotate_ch = fasta_top_ch
            .map    { gid, fa -> tuple(gid, 'top', fa) }
            .combine(top_species_ch)
            .combine(tree_ch)

        def bottom_annotate_ch = fasta_bottom_ch
            .map    { gid, fa -> tuple(gid, 'bottom', fa) }
            .combine(bottom_species_ch)
            .combine(tree_ch)

        def annotate_input_ch = top_annotate_ch.mix(bottom_annotate_ch)
            // → (gene_id, direction, fasta, fg_species_file, tree)

        def annotBatchSize = (params.fade_batch_size ?: 1) as int
        def annotated_ch
        def filtered_fasta_ch

        if (annotBatchSize > 1) {
            // ── Batched annotation mode ──────────────────────────────────────
            def annotate_branched = annotate_input_ch.branch {
                top:    it[1] == 'top'
                bottom: it[1] == 'bottom'
            }

            def make_annotate_batches = { branch_ch, dir ->
                def batchCounter = 0
                branch_ch
                    .toSortedList({ a, b -> a[0] <=> b[0] })
                    .flatMap()
                    .collate(annotBatchSize)
                    .map { batch ->
                        def batchID = sprintf('annotate_batch_%s_%05d', dir, ++batchCounter)
                        def validRows = batch.findAll { row ->
                            row[2]?.name != 'NO_FILE' && row[3]?.name != 'NO_FILE' && row[4]?.name != 'NO_FILE'
                        }
                        if (!validRows) return null
                        def manifestText = createBatchManifestText(
                            validRows.collect { row -> "${row[0]}\t${row[1]}\t${row[2].name}\t${row[4].name}\t${row[3].name}" }
                        )
                        // Tree and FG list are shared within a direction batch;
                        // deduplicate to avoid Nextflow stageAs filename collisions.
                        def uniqTrees = validRows.collect { row -> row[4] }.unique { it.name }
                        def uniqSpecies = validRows.collect { row -> row[3] }.unique { it.name }
                        tuple(batchID, validRows.size(), manifestText,
                              validRows.collect { row -> row[2] },   // fastas
                              uniqTrees,                          // trees
                              uniqSpecies)                        // species files
                    }
                    .filter { it != null }
            }

            def batches_ch = make_annotate_batches(annotate_branched.top,    'top')
                .mix(make_annotate_batches(annotate_branched.bottom, 'bottom'))

            def batched_out = ANNOTATE_TREE_FG_BATCHED(batches_ch)
            
            // Reconstruct (gene_id, direction, file) tuples from output filenames
            annotated_ch = batched_out.annotated_trees
                .flatten()
                .map { f ->
                    def m = (f.name =~ /^(.+)_(top|bottom)_fg\.nwk$/)
                    tuple(m[0][1], m[0][2], f)
                }

            filtered_fasta_ch = batched_out.filtered_fastas
                .flatten()
                .map { f ->
                    def m = (f.name =~ /^(.+)_(top|bottom)\.fa$/)
                    tuple(m[0][1], m[0][2], f)
                }
        } else {
            // ── Single-gene annotation mode ──────────────────────────────────
            def annotate_result = ANNOTATE_TREE_FG(annotate_input_ch)
            annotated_ch = annotate_result.annotated_tree
            filtered_fasta_ch = annotate_result.filtered_fasta
        }

        // ── Join filtered FASTA + annotated tree, then run FADE ─────────────
        //  filtered_fasta_ch : (gene_id, direction, filtered_fasta)
        //  annotated_ch      : (gene_id, direction, annotated_tree)
        //  join by [0,1]     → (gene_id, direction, filtered_fasta, annotated_tree)
        def fade_input_ch = filtered_fasta_ch.join(annotated_ch, by: [0, 1])

        def fadeBatchSize = (params.fade_batch_size ?: 1) as int
        if (fadeBatchSize > 1) {
            // ── Batched mode ──────────────────────────────────────────────
            def fade_branched = fade_input_ch.branch {
                top:    it[1] == 'top'
                bottom: it[1] == 'bottom'
            }

            def make_fade_batches = { branch_ch, dir ->
                def batchCounter = 0
                branch_ch
                    .toSortedList({ a, b -> a[0] <=> b[0] })
                    .flatMap()
                    .collate(fadeBatchSize)
                    .map { batch ->
                        def batchID = sprintf('fade_batch_%s_%05d', dir, ++batchCounter)
                        def manifestText = createBatchManifestText(
                            batch.collect { row -> "${row[0]}\t${row[2].name}\t${row[3].name}" }
                        )
                        tuple(batchID, dir, batch.size(), manifestText,
                              batch.collect { row -> row[2] },   // filtered fastas
                              batch.collect { row -> row[3] })   // annotated trees
                    }
            }

            def batches_ch = make_fade_batches(fade_branched.top,    'top')
                .mix(make_fade_batches(fade_branched.bottom, 'bottom'))

            def batched_out = FADE_BATCHED(batches_ch, lg_dat_ch).fade_json
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
        def branched = fade_results_ch.branch {
            top:    it[1] == 'top'
            bottom: it[1] == 'bottom'
        }
        def top_jsons    = branched.top.map    { it[2] }.collect().ifEmpty([])
        def bottom_jsons = branched.bottom.map { it[2] }.collect().ifEmpty([])

        def fg_top_ch    = top_species_ch.ifEmpty    { file('NO_FG_LIST') }
        def fg_bottom_ch = bottom_species_ch.ifEmpty { file('NO_FG_LIST') }

        fade_report_top    = FADE_REPORT_TOP(Channel.value('top'),       top_jsons,    fg_top_ch   )
        fade_report_bottom = FADE_REPORT_BOTTOM(Channel.value('bottom'), bottom_jsons, fg_bottom_ch)

        // ── Position-level site CSV (always runs — posenrich's FADE evidence
        // layer, independent of --enrichment/--string) ───────────────────────
        def fade_sites_top    = FADE_JSON_TO_CSV_TOP(Channel.value('top'),       top_jsons)
        def fade_sites_bottom = FADE_JSON_TO_CSV_BOTTOM(Channel.value('bottom'), bottom_jsons)

        // ── Gene list extraction + FCS/STRING per direction ──────────────────
        if (params.enrichment || params.string) {
            def run_fade_enrich = { direction, summary_tsv_ch, gene_lists_proc ->
                def lists_out    = gene_lists_proc(Channel.value(direction), summary_tsv_ch)
                def subpath      = "selection/fade/${direction}"
                def bg_ch        = lists_out.gene_lists.flatten().filter { it.name == 'background.txt' }.collect()
                def interest_ch  = lists_out.gene_lists.flatten().filter { it.name != 'background.txt' }.collect()
                // FADE no longer runs its own FCS ranking (see fcs.nf/scoring_enrichment.nf —
                // its max-of-many-sites statistic proved unreliable for gene-set enrichment;
                // FADE now contributes as cross-module corroboration + a posenrich position
                // group instead). STRING (set-based, not rank-based) is unaffected.
                if (params.string) {
                    FADE_STRING_REPORT(
                        Channel.value(subpath), bg_ch, interest_ch,
                        Channel.value("15.STRING_fade_${direction}_${params.traitname ?: 'trait'}")
                    )
                }
            }

            run_fade_enrich('top',    fade_report_top.summary_tsv,    { d, s -> FADE_GENE_LISTS_TOP(d, s) })
            run_fade_enrich('bottom', fade_report_bottom.summary_tsv, { d, s -> FADE_GENE_LISTS_BOTTOM(d, s) })
        }

    emit:
        report_top     = fade_report_top.report
        report_bottom  = fade_report_bottom.report
        summary_top    = fade_report_top.summary_tsv
        summary_bottom = fade_report_bottom.summary_tsv
        site_tsv_top    = fade_report_top.site_tsv
        site_tsv_bottom = fade_report_bottom.site_tsv
        json_results   = fade_results_ch
        sites_csv_top    = fade_sites_top.sites_csv
        sites_csv_bottom = fade_sites_bottom.sites_csv
}
