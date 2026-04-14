#!/usr/bin/env nextflow

/*
#
#  ███╗   ███╗ ██████╗ ██╗     ███████╗██████╗  █████╗ ████████╗███████╗
#  ████╗ ████║██╔═══██╗██║     ██╔════╝██╔══██╗██╔══██╗╚══██╔══╝██╔════╝
#  ██╔████╔██║██║   ██║██║     █████╗  ██████╔╝███████║   ██║   █████╗
#  ██║╚██╔╝██║██║   ██║██║     ██╔══╝  ██╔══██╗██╔══██║   ██║   ██╔══╝
#  ██║ ╚═╝ ██║╚██████╔╝███████╗███████╗██║  ██║██║  ██║   ██║   ███████╗
#  ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   ╚══════╝
#
# PHYLOPHERE: MoleRate Workflow
#
# Tests foreground branches (defined by trait_stats.csv extremes) for
# differential evolutionary rate using HyPhy MoleRate
# (hyphy-analyses/MoleRate/MoleRate.bf).
#
# Two directions always run in parallel:
#   top    -> global_label == high_extreme species as foreground
#   bottom -> global_label == low_extreme species as foreground
#
# Alignment preparation (PHYLIP→FASTA conversion, tree-filtering, and species
# extraction) is handled upstream by SELECTION_PREP (called once in main.nf),
# which shares the results with FADE to avoid redundant processing.
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/molerate.nf
*/

include { EXTRACT_FG_BRANCHES             } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { MOLERATE_RUN; MOLERATE_BATCHED  } from "${baseDir}/subworkflows/MOLERATE/molerate_run.nf"
include { MOLERATE_REPORT as MOLERATE_REPORT_TOP    } from "${baseDir}/subworkflows/MOLERATE/molerate_report.nf"
include { MOLERATE_REPORT as MOLERATE_REPORT_BOTTOM } from "${baseDir}/subworkflows/MOLERATE/molerate_report.nf"


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

workflow MOLERATE {

    take:
        // Pre-processed alignment channels from SELECTION_PREP (called once in main.nf).
        // Each channel emits (gene_id, filtered_fasta) tuples — no direction label yet;
        // direction is introduced inside this workflow when building the MoleRate inputs.
        fasta_top_ch      // (gene_id, fasta) for top-direction genes
        fasta_bottom_ch   // (gene_id, fasta) for bottom-direction genes
        // Shared reference files (value channels — can be broadcast to every gene)
        tree_ch           // path to phenotype-pruned species tree
        top_species_ch    // path to top_species.txt
        bottom_species_ch // path to bottom_species.txt

    main:

        // ── LG model dat file ────────────────────────────────────────────────
        def lg_dat_ch = Channel.value(file(params.lg_dat_path))

        // ── Extract FG branch lists per direction ────────────────────────────
        // EXTRACT_FG_BRANCHES converts a species-file into a plain taxon-name
        // list consumed by MoleRate's --branches argument.
        def fg_inputs_ch = top_species_ch
            .map { f -> tuple('top', f) }
            .mix(bottom_species_ch.map { f -> tuple('bottom', f) })

        fg_branches_ch = EXTRACT_FG_BRANCHES(fg_inputs_ch).fg_list
        // -> (direction, fg_list_file)

        // ── Build per-direction MoleRate input channel ───────────────────────
        // Introduce direction labels and join with fg_branches and tree.
        //
        // fasta_top_ch    : (gene_id, fasta)
        // fasta_bottom_ch : (gene_id, fasta)
        // fg_branches_ch  : (direction, fg_list_file)
        // tree_ch         : path (value — broadcasts to every gene)
        //
        // Strategy: tag each fasta with its direction, then combine with
        // fg_branches_ch on direction (broadcast 2 fg_list files across all genes),
        // then combine with tree_ch (broadcast single tree).

        def top_keyed_ch    = fasta_top_ch.map    { gid, fa -> tuple('top',    gid, fa) }
        def bottom_keyed_ch = fasta_bottom_ch.map { gid, fa -> tuple('bottom', gid, fa) }
        def fasta_keyed_ch  = top_keyed_ch.mix(bottom_keyed_ch)
        // -> (direction, gene_id, fasta)

        molerate_input_ch = fasta_keyed_ch
            .combine(fg_branches_ch, by: 0)
            // -> (direction, gene_id, fasta, fg_list)
            .combine(tree_ch)
            // -> (direction, gene_id, fasta, fg_list, tree)
            .map { dir, gid, fa, fg_list, tree -> tuple(gid, dir, fa, tree, fg_list) }
            // -> (gene_id, direction, fasta, tree, fg_list)  [matches MOLERATE_RUN input order]

        def molerateBatchSize = (params.molerate_batch_size ?: 1) as int
        if (molerateBatchSize > 1) {
            // ── Batched mode ──────────────────────────────────────────────
            // tree and fg_list are shared per direction — take from batch[0].
            def molerate_branched = molerate_input_ch.branch {
                top:    it[1] == 'top'
                bottom: it[1] == 'bottom'
            }

            def make_molerate_batches = { branch_ch, dir ->
                branch_ch
                    .collate(molerateBatchSize)
                    .map { batch ->
                        def batchID = "molerate_batch_${dir}_${java.util.UUID.randomUUID().toString().replace('-', '').take(12)}"
                        def validRows = batch.findAll { row ->
                            row[2]?.name != 'NO_FILE' && row[3]?.name != 'NO_FILE' && row[4]?.name != 'NO_FILE'
                        }
                        if (!validRows) return null
                        def manifestText = createBatchManifestText(
                            validRows.collect { row -> "${row[0]}\t${row[2].name}" }
                        )
                        def sharedTree   = validRows[0][3]
                        def sharedFgList = validRows[0][4]
                        tuple(batchID, dir, validRows.size(), manifestText,
                              validRows.collect { row -> row[2] },   // fastas
                              sharedTree, sharedFgList)
                    }
                    .filter { it != null }
            }

            def batches_ch = make_molerate_batches(molerate_branched.top,    'top')
                .mix(make_molerate_batches(molerate_branched.bottom, 'bottom'))

            def batched_out = MOLERATE_BATCHED(batches_ch, lg_dat_ch).molerate_json
            molerate_results_ch = batched_out
                .transpose()
                .map { dir, f ->
                    def m = (f.name =~ /^(.+)\.(top|bottom)\.molerate\.json$/)
                    tuple(m[0][1], m[0][2], f)
                }
        } else {
            // ── Single-gene mode (default) ────────────────────────────────
            molerate_results_ch = MOLERATE_RUN(molerate_input_ch, lg_dat_ch).molerate_json
        }

        // ── Reports per direction ────────────────────────────────────────────
        def branched = molerate_results_ch.branch {
            top:    it[1] == 'top'
            bottom: it[1] == 'bottom'
        }
        def top_jsons    = branched.top.map    { it[2] }.collect().ifEmpty([])
        def bottom_jsons = branched.bottom.map { it[2] }.collect().ifEmpty([])

        def fg_top_ch    = fg_branches_ch.filter { it[0] == 'top'    }.map { it[1] }
            .ifEmpty { file('NO_FG_LIST') }
        def fg_bottom_ch = fg_branches_ch.filter { it[0] == 'bottom' }.map { it[1] }
            .ifEmpty { file('NO_FG_LIST') }

        def rpt_top    = MOLERATE_REPORT_TOP(Channel.value('top'),       top_jsons,    fg_top_ch   )
        def rpt_bottom = MOLERATE_REPORT_BOTTOM(Channel.value('bottom'), bottom_jsons, fg_bottom_ch)

    emit:
        report_top     = rpt_top.report
        report_bottom  = rpt_bottom.report
        summary_top    = rpt_top.summary_tsv
        summary_bottom = rpt_bottom.summary_tsv
        json_results   = molerate_results_ch
}
