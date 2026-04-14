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
# Tests foreground branches (defined by the caastools traitfile) for
# differential evolutionary rate using HyPhy MoleRate
# (hyphy-analyses/MoleRate/MoleRate.bf).
#
# Two directions always run in parallel:
#   top    -> contrast_group == 1 species as foreground
#   bottom -> contrast_group == 0 species as foreground
#
# Run modes (params.molerate_mode):
#   gene_set  -> only genes from CT accumulation / postproc significant lists
#   all       -> every PHYLIP alignment file in the alignment directory
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: workflows/molerate.nf
*/

include { COLLECT_GENE_SETS     } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { PHYLIP_TO_FASTA       } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { EXTRACT_FG_BRANCHES   } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { MOLERATE_RUN          } from "${baseDir}/subworkflows/MOLERATE/molerate_run.nf"
include { MOLERATE_REPORT       } from "${baseDir}/subworkflows/MOLERATE/molerate_report.nf"


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
        def gid = f.name.replaceAll(/\.(phy|phylip|aln)$/, '')
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
        traitfile_input   // Channel<path> or null -> falls back to params.caas_config
        tree_input        // Channel<path> or null -> falls back to params.tree
        // Optional upstream channels for inline gene-set piping.
        // Pass Channel.empty() when running standalone (params-based).
        acc_top_ch        // accumulation CSV for TOP  (*_top_aggregated_results.csv)
        acc_bottom_ch     // accumulation CSV for BOTTOM
        pp_top_ch         // postproc TXT for TOP  (*_change_side_top_significant.txt)
        pp_bottom_ch      // postproc TXT for BOTTOM

    main:

        // 1. Resolve shared inputs ----------------------------------------
        def traitfile_ch = (traitfile_input ?: Channel.empty())
            .ifEmpty {
                assert params.caas_config : \
                    "MoleRate requires a CONTRAST_SELECTION traitfile or --caas_config"
                def f = file(params.caas_config)
                assert f.exists() : "MoleRate: traitfile not found: ${params.caas_config}"
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

        // 2. Extract FG branch lists per direction (one process call each) ---
        //    EXTRACT_FG_BRANCHES input: (direction, traitfile)
        def fg_inputs_ch = Channel.from(['top', 'bottom'])
            .combine(traitfile_ch)
        // -> (direction, traitfile)

        fg_branches_ch = EXTRACT_FG_BRANCHES(fg_inputs_ch).fg_branches
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
            def resolved_acc_top = (acc_top_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_accumulation_top    ?: 'NO_FILE') }

            def resolved_acc_bottom = (acc_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_accumulation_bottom ?: 'NO_FILE') }

            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_postproc_top        ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.molerate_postproc_bottom     ?: 'NO_FILE') }

            def sets = COLLECT_GENE_SETS(
                resolved_acc_top,
                resolved_acc_bottom,
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
        //  Strategy: cross-join fasta_ch with fg_branches_ch on direction, then
        //  combine with tree_ch (single value).
        //
        //  fasta_ch.map -> (direction, gene_id, fasta)          [keyed by direction]
        //  join fg_branches_ch -> (direction, gene_id, fasta, fg_list)
        //  combine tree_ch -> (direction, gene_id, fasta, fg_list, tree)

        def fasta_keyed_ch = fasta_ch.map { gid, dir, fa -> tuple(dir, gid, fa) }

        molerate_input_ch = fasta_keyed_ch
            .join(fg_branches_ch, by: 0)
            // -> (direction, gene_id, fasta, fg_list)
            .combine(tree_ch)
            // -> (direction, gene_id, fasta, fg_list, tree)
            .map { dir, gid, fa, fg_list, tree -> tuple(gid, dir, fa, fg_list, tree) }
            // -> (gene_id, direction, fasta, fg_list, tree)

        molerate_results_ch = MOLERATE_RUN(molerate_input_ch, lg_dat_ch).molerate_json

        // 6. Reports per direction -----------------------------------------
        def top_jsons = molerate_results_ch
            .filter { gid, dir, json -> dir == 'top' }
            .map    { gid, dir, json -> json }
            .collect()

        def bottom_jsons = molerate_results_ch
            .filter { gid, dir, json -> dir == 'bottom' }
            .map    { gid, dir, json -> json }
            .collect()

        def rpt_top    = MOLERATE_REPORT(Channel.value('top'),    top_jsons   )
        def rpt_bottom = MOLERATE_REPORT(Channel.value('bottom'), bottom_jsons)

    emit:
        report_top     = rpt_top.report
        report_bottom  = rpt_bottom.report
        summary_top    = rpt_top.summary_tsv
        summary_bottom = rpt_bottom.summary_tsv
        json_results   = molerate_results_ch
}
