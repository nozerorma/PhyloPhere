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
# Run modes (params.molerate_mode):
#   gene_set  -> only genes from CT postproc significant lists
#   all       -> every PHYLIP alignment file in the alignment directory
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/molerate.nf
*/

include { COLLECT_GENE_SETS       } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { EXTRACT_EXTREME_SPECIES } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { PHYLIP_TO_FASTA         } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { EXTRACT_FG_BRANCHES     } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { MOLERATE_RUN; MOLERATE_BATCHED } from "${baseDir}/subworkflows/MOLERATE/molerate_run.nf"
include { MOLERATE_REPORT as MOLERATE_REPORT_TOP    } from "${baseDir}/subworkflows/MOLERATE/molerate_report.nf"
include { MOLERATE_REPORT as MOLERATE_REPORT_BOTTOM } from "${baseDir}/subworkflows/MOLERATE/molerate_report.nf"


// ---------------------------------------------------------------------------
// Render a list of TSV row strings into a heredoc-safe manifest block.
def createBatchManifestText = { List<String> rows ->
    rows
        .collect { row -> row.replaceFirst(/^\s+/, '') }
        .join(System.lineSeparator()) + System.lineSeparator()
}

// ---------------------------------------------------------------------------
// Private helper: build a List of [gene_id, direction, phylip_File] tuples.
// Runs at channel-operator evaluation time (inside flatMap).
// ---------------------------------------------------------------------------
def ali_tuples_from_dir(String ali_dir, String direction, Set wanted) {
    def ali_path = file(ali_dir)
    if (!ali_path.isDirectory()){
        log.warn "MoleRate: alignment directory not found: ${ali_dir} -- skipping ${direction}"
        return []
    }
    def all_files = ali_path.listFiles()?.findAll { f ->
        f.isFile() && (
            f.name.endsWith('.phy')    ||
            f.name.endsWith('.phylip') ||
            f.name.endsWith('.aln')    ||
            !f.name.contains('.')
        )
    } ?: []

    def tuples = all_files.collect { f ->
        def gid = f.name.tokenize('.')[0]
        (wanted == null || wanted.contains(gid)) ? [gid, direction, f] : null
    }.findAll { it != null }

    if (tuples.isEmpty()) {
        log.warn "MoleRate: no alignment files matched for direction '${direction}' in ${ali_dir}"
    } else {
        log.info "MoleRate: ${tuples.size()} genes queued for '${direction}' direction"
    }
    return tuples
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow MOLERATE {

    take:
        stats_file_input  // Channel<path> or null -> falls back to canonical trait_stats.csv under outdir
        tree_input        // Channel<path> or null -> falls back to params.tree
        // Optional upstream channels for inline gene-set piping.
        // Pass Channel.empty() when running standalone (params-based).
        pp_top_ch         // postproc TXT for TOP  (all_top.txt)
        pp_bottom_ch      // postproc TXT for BOTTOM (all_bottom.txt)

    main:

        // 1. Resolve shared inputs ----------------------------------------
        def stats_file_ch = (stats_file_input ?: Channel.empty())
            .ifEmpty {
                def fallback = "${params.outdir}/data_exploration/1.Data-exploration/1.Species_distribution/trait_stats.csv"
                def f = file(fallback)
                assert f.exists() : "MoleRate: trait_stats.csv not found at ${fallback}"
                f
            }

        def tree_ch = (tree_input ?: Channel.empty())
            .ifEmpty {
                assert params.tree : "MoleRate requires --tree"
                def f = file(params.tree)
                assert f.exists() : "MoleRate: tree not found: ${params.tree}"
                f
            }

        def ali_dir = (params.molerate_alignment && params.molerate_alignment != '') \
            ? params.molerate_alignment : params.alignment
        assert ali_dir : \
            "MoleRate: no alignment directory specified (--molerate_alignment or --alignment)"

        // LG model dat file
        def lg_dat_ch = Channel.value(file(params.lg_dat_path))
        def species_lists = EXTRACT_EXTREME_SPECIES(stats_file_ch)
        def fg_inputs_ch = species_lists.top_species
            .map { f -> tuple('top', f) }
            .mix(species_lists.bottom_species.map { f -> tuple('bottom', f) })

        // 2. Extract FG branch lists per direction (one process call each) ---
        //    EXTRACT_FG_BRANCHES input: (direction, fg_species_file)

        fg_branches_ch = EXTRACT_FG_BRANCHES(fg_inputs_ch).fg_list
        // -> (direction, fg_list_file)

        // 3. Build per-direction alignment channels -----------------------
        def all_ali_ch

        if (params.molerate_mode == 'all') {
            def top_tuples    = ali_tuples_from_dir(ali_dir, 'top',    null)
            def bottom_tuples = ali_tuples_from_dir(ali_dir, 'bottom', null)
            all_ali_ch = Channel.fromList(top_tuples + bottom_tuples)

        } else {
            // gene_set mode -----------------------------------------------
            // Prefer upstream piped channels; fall back to --molerate_* path params.
            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_postproc_top        ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_postproc_bottom     ?: 'NO_FILE') }

            def sets = COLLECT_GENE_SETS(
                resolved_pp_top,
                resolved_pp_bottom
            )

            def top_ali_ch = sets.gene_set_top.flatMap { gsf ->
                def wanted = gsf.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#')}.toSet()
                ali_tuples_from_dir(ali_dir, 'top', wanted)
            }

            def bottom_ali_ch = sets.gene_set_bottom.flatMap { gsf ->
                def wanted = gsf.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#')}.toSet()
                ali_tuples_from_dir(ali_dir, 'bottom', wanted)
            }

            all_ali_ch = top_ali_ch.mix(bottom_ali_ch)
        }

        // 4. PHYLIP -> FASTA conversion ------------------------------------
        fasta_ch = PHYLIP_TO_FASTA(all_ali_ch).fasta
        // -> (gene_id, direction, fasta)

        // 5. Join fasta + fg_branches_ch by direction, broadcast tree -----
        //  fasta_ch      : (gene_id, direction, fasta)
        //  fg_branches_ch: (direction, fg_list)
        //  For each (gene_id, direction) we need: fasta, fg_list, tree.
        //
        //  Strategy: combine (broadcast) fasta_ch with fg_branches_ch on direction, then
        //  combine with tree_ch (single value).
        //
        //  fasta_ch.map -> (direction, gene_id, fasta)          [keyed by direction]
        //  combine fg_branches_ch -> (direction, gene_id, fasta, fg_list)
        //                           [broadcasts the 2 fg_list files across all genes]
        //  combine tree_ch -> (direction, gene_id, fasta, fg_list, tree)

        def fasta_keyed_ch = fasta_ch.map { gid, dir, fa -> tuple(dir, gid, fa) }

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
            // molerate_input_ch: (gene_id, direction, fasta, tree, fg_list)
            // tree and fg_list are shared per direction — take from batch[0].
            def molerate_branched = molerate_input_ch.branch {
                top:    it[1] == 'top'
                bottom: it[1] == 'bottom'
            }
            def molerateBatchCounterTop    = 0
            def molerateBatchCounterBottom = 0

            def make_molerate_batches = { branch_ch, dir, counter ->
                branch_ch
                    .collate(molerateBatchSize)
                    .map { batch ->
                        def batchID = sprintf("molerate_batch_${dir}_%05d", ++counter)
                        def manifestText = createBatchManifestText(
                            batch.collect { row -> "${row[0]}\t${row[2].name}" }
                        )
                        def sharedTree   = batch[0][3]
                        def sharedFgList = batch[0][4]
                        tuple(batchID, dir, batch.size(), manifestText,
                              batch.collect { row -> row[2] },   // fastas
                              sharedTree, sharedFgList)
                    }
            }

            def batches_ch = make_molerate_batches(molerate_branched.top,    'top',    molerateBatchCounterTop)
                .mix(make_molerate_batches(molerate_branched.bottom, 'bottom', molerateBatchCounterBottom))

            def batched_out = MOLERATE_BATCHED(batches_ch, lg_dat_ch).molerate_json
            // molerate_json emits (direction, path) per JSON; flatten to (gene_id, direction, path)
            molerate_results_ch = batched_out
                .transpose()
                .map { dir, f ->
                    def gene_id = f.name.replace(".${dir}.molerate.json", "")
                    tuple(gene_id, dir, f)
                }
        } else {
            // ── Single-gene mode (default) ────────────────────────────────
            molerate_results_ch = MOLERATE_RUN(molerate_input_ch, lg_dat_ch).molerate_json
        }

        // 6. Reports per direction -----------------------------------------
        //    Use .branch() to split once and guarantee exclusive routing,
        //    avoiding any ambiguity from applying two independent .filter()
        //    subscriptions on the same queue channel.
        def branched = molerate_results_ch.branch {
            top:    it[1] == 'top'
            bottom: it[1] == 'bottom'
        }
        def top_jsons    = branched.top.map    { it[2] }.collect().ifEmpty([])
        def bottom_jsons = branched.bottom.map { it[2] }.collect().ifEmpty([])

        def rpt_top    = MOLERATE_REPORT_TOP(Channel.value('top'),    top_jsons   )
        def rpt_bottom = MOLERATE_REPORT_BOTTOM(Channel.value('bottom'), bottom_jsons)

    emit:
        report_top     = rpt_top.report
        report_bottom  = rpt_bottom.report
        summary_top    = rpt_top.summary_tsv
        summary_bottom = rpt_bottom.summary_tsv
        json_results   = molerate_results_ch
}
