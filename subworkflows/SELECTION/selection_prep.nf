#!/usr/bin/env nextflow

/*
 * SELECTION / shared alignment preprocessing
 *
 * SELECTION_PREP runs the alignment preparation pipeline ONCE and shares the
 * results between FADE and MOLERATE, eliminating duplicate work:
 *
 *   1. EXTRACT_EXTREME_SPECIES  – top/bottom species lists from trait_stats.csv
 *   2. [Alignment channel]      – build gene list (all / gene_set / toy_mode)
 *   3. PREP_ALIGNMENTS_BATCHED  – convert PHYLIP→FASTA + filter to tree taxa
 *                                  (or PREP_ALIGNMENTS for batch_size=1)
 *
 * Outputs:
 *   filtered_fasta_top_ch    (gene_id, fasta) for top-direction genes
 *   filtered_fasta_bottom_ch (gene_id, fasta) for bottom-direction genes
 *   top_species              path to top_species.txt   (value channel)
 *   bottom_species           path to bottom_species.txt (value channel)
 *   tree_ch                  path to pruned tree        (value channel)
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 * File: subworkflows/SELECTION/selection_prep.nf
 */

include { EXTRACT_EXTREME_SPECIES } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"
include { COLLECT_GENE_SETS       } from "${baseDir}/subworkflows/SELECTION/selection_utils.nf"


// ─── Helpers ─────────────────────────────────────────────────────────────────

/**
 * Render a list of TSV row strings into a heredoc-safe manifest block.
 */
def createBatchManifestText = { List<String> rows ->
    rows
        .collect { row -> row.replaceFirst(/^\s+/, '') }
        .join(System.lineSeparator()) + System.lineSeparator()
}

/**
 * Build a List of [gene_id, alignment_File] tuples synchronously (no direction).
 * Returns a plain Groovy List (not a Channel) so it is safe to call inside a
 * channel operator's closure (e.g. flatMap).
 *
 *   wanted  — Set<String> of gene IDs to keep, or null / empty to keep all.
 */
def ali_tuples_from_dir(String ali_dir, Set wanted) {
    def ali_path = file(ali_dir)
    if (!ali_path.isDirectory()) {
        log.warn "SELECTION_PREP: alignment directory not found: ${ali_dir}"
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
        (wanted == null || wanted.contains(gid)) ? [gid, f] : null
    }.findAll { it != null }

    if (tuples.isEmpty())
        log.warn "SELECTION_PREP: no alignment files matched in ${ali_dir}"
    else
        log.info  "SELECTION_PREP: ${tuples.size()} genes queued for preprocessing"
    return tuples
}


// ─── Processes ───────────────────────────────────────────────────────────────

/**
 * PREP_ALIGNMENTS — single-gene fallback (batch_size = 1).
 * Converts alignment to FASTA (if needed) then restricts sequences to the
 * taxa present in the phenotype-pruned tree.
 * Output file: ${gene_id}.filtered.fa  (no direction — filtering is direction-agnostic)
 */
process PREP_ALIGNMENTS {
    tag "${gene_id}"

    input:
    tuple val(gene_id), path(alignment_file), path(tree)

    output:
    tuple val(gene_id), path("${gene_id}.filtered.fa"), emit: filtered_fasta, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    def pyBin = (params.use_singularity || params.use_apptainer)
        ? '/usr/local/bin/_entrypoint.sh python'
        : 'python'
    """
    # Convert to FASTA (or copy if already FASTA)
    _lower=\$(echo "${alignment_file}" | tr '[:upper:]' '[:lower:]')
    if [[ "\$_lower" == *.fa || "\$_lower" == *.fasta ]]; then
        cp "${alignment_file}" "${gene_id}_tmp.fa"
    else
        ${pyBin} ${local_dir}/phylip_to_fasta.py "${alignment_file}" "${gene_id}_tmp.fa"
    fi

    # Filter to tree taxa
    ${pyBin} ${local_dir}/filter_fasta_to_tree.py \\
        --tree "${tree}" \\
        --fasta "${gene_id}_tmp.fa" \\
        --output "${gene_id}.filtered.fa"

    [ -s "${gene_id}.filtered.fa" ] || rm -f "${gene_id}.filtered.fa"
    """
}


/**
 * PREP_ALIGNMENTS_BATCHED — batched conversion + tree-filtering.
 * Processes multiple genes in a single Nextflow task, reducing SLURM scheduling
 * overhead. Up to params.selection_prep_batch_workers genes run concurrently
 * within each task using bash job control.
 *
 * Manifest format (tab-separated, one gene per line):
 *   gene_id <TAB> alignment_filename
 *
 * Alignment files are staged under:   alignments/<alignment_filename>
 * Output files:                        <gene_id>.filtered.fa
 */
process PREP_ALIGNMENTS_BATCHED {
    tag "${batchID}"

    input:
    tuple val(batchID), val(batchSize), val(manifestText),
          path(alignments, stageAs: 'alignments/*')
    path tree

    output:
    tuple val(batchID), path("*.filtered.fa"), emit: filtered_fastas, optional: true

    script:
    def local_dir  = "${baseDir}/subworkflows/SELECTION/local"
    def script_dir = "${baseDir}/subworkflows/SELECTION/local/scripts"
    def workers    = params.selection_prep_batch_workers ?: 4
    def runnerMode = (params.use_singularity || params.use_apptainer) ? "container" : "local"
    """
    cat <<'MANIFEST_EOF' > ${batchID}.manifest.tsv
${manifestText}MANIFEST_EOF

    bash ${script_dir}/run_prep_alignments_batch.sh \\
        --batch-id    "${batchID}" \\
        --manifest    "${batchID}.manifest.tsv" \\
        --workers     "${workers}" \\
        --runner-mode "${runnerMode}" \\
        --tree        "${tree}" \\
        --local-dir   "${local_dir}"
    """
}


// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow SELECTION_PREP {

    take:
        stats_file_input   // Channel<path> or null → falls back to canonical outdir path
        tree_input         // Channel<path> or null → falls back to params.tree
        pp_top_ch          // postproc TXT for TOP  (all_top.txt)  — gene_set mode
        pp_bottom_ch       // postproc TXT for BOTTOM (all_bottom.txt) — gene_set mode
        ct_discovery_input // Optional CT discovery table for toy_mode gene reuse

    main:

        // ── Resolve shared inputs ────────────────────────────────────────────
        // Convert to value channels so they broadcast correctly to every gene.
        tree_val = (tree_input ?: Channel.empty())
            .ifEmpty {
                assert params.tree : "SELECTION_PREP requires --tree"
                // Tree is a user-provided parameter and must exist.
                def f = file(params.tree)
                assert f.exists() : "SELECTION_PREP: tree not found: ${params.tree}"
                f
            }
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files -> files[0] }

        def stats_file_ch = (stats_file_input ?: Channel.empty())
            .ifEmpty {
                def fallback = "${params.outdir}/data_exploration/1.Data-exploration/1.Species_distribution/trait_stats.csv"
                def f = file(fallback)
                // Don't assert existence here — file is created by upstream DATASET_EXPLORATION process
                f
            }
            .collect()
            .filter { files -> files && files.size() > 0 }
            .map { files -> files[0] }

        def ali_dir = (params.fade_alignment ?: params.molerate_alignment ?: params.alignment) ?:
            error("SELECTION_PREP: no alignment directory specified (--alignment, --fade_alignment, or --molerate_alignment)")

        // Prefer the FADE mode if set; fall back to MOLERATE mode, then 'gene_set'.
        def mode = (params.fade_mode ?: params.molerate_mode ?: 'gene_set')

        // ── Extract foreground species (once, shared by FADE and MOLERATE) ───
        def species_lists      = EXTRACT_EXTREME_SPECIES(stats_file_ch)
        top_species_val        = species_lists.top_species.collect().map { it[0] }
        bottom_species_val     = species_lists.bottom_species.collect().map { it[0] }

        // ── Build alignment channel (gene_id, alignment_file) ────────────────
        // No direction at this stage — conversion and tree-filtering are
        // direction-agnostic; direction is introduced later inside FADE/MOLERATE.

        def ali_ch
        // Gene-set files are only needed in gene_set mode; initialise as empty.
        def gs_top_val    = Channel.value(file('NO_FILE'))
        def gs_bottom_val = Channel.value(file('NO_FILE'))

        if (mode == 'all') {
            if (params.toy_mode) {
                def resolved_ct = (ct_discovery_input ?: Channel.empty())
                    .ifEmpty { file('NO_FILE') }

                ali_ch = resolved_ct.flatMap { discovery_file ->
                    def toy_genes = null
                    if (discovery_file.name != 'NO_FILE' && discovery_file.exists()) {
                        toy_genes = discovery_file
                            .readLines()
                            .drop(1)
                            .collect { it.split('\t')[0].trim() }
                            .findAll { it }
                            .toSet()
                        log.info "[toy_mode] SELECTION_PREP: reusing ${toy_genes.size()} genes from CT discovery"
                    } else {
                        def n = (params.toy_n ?: 50) as int
                        def all_tuples = ali_tuples_from_dir(ali_dir, null)
                        Collections.shuffle(all_tuples)
                        toy_genes = all_tuples.take(n).collect { it[0] }.toSet()
                        log.info "[toy_mode] SELECTION_PREP: using ${toy_genes.size()} randomly sampled genes"
                    }
                    ali_tuples_from_dir(ali_dir, toy_genes)
                }
            } else {
                ali_ch = Channel.fromList(ali_tuples_from_dir(ali_dir, null))
            }

        } else {
            // gene_set mode ───────────────────────────────────────────────────
            // Prefer upstream piped channels; fall back to --fade_* / --molerate_* params.
            def resolved_pp_top = (pp_top_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_top ?: params.molerate_postproc_top ?: 'NO_FILE') }

            def resolved_pp_bottom = (pp_bottom_ch ?: Channel.empty())
                .ifEmpty { file(params.fade_postproc_bottom ?: params.molerate_postproc_bottom ?: 'NO_FILE') }

            def gene_sets = COLLECT_GENE_SETS(resolved_pp_top, resolved_pp_bottom)

            // Keep references for direction-routing below.
            gs_top_val    = gene_sets.gene_set_top.collect().map { it[0] }
            gs_bottom_val = gene_sets.gene_set_bottom.collect().map { it[0] }

            // Preprocess the UNION of top+bottom gene sets so each gene is
            // converted and tree-filtered exactly once.
            ali_ch = gene_sets.gene_set_top.combine(gene_sets.gene_set_bottom).flatMap { gsf_top, gsf_bottom ->
                def top_genes = gsf_top.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#') }.toSet()
                def bottom_genes = gsf_bottom.readLines()
                    .collect { it.trim() }.findAll { it && !it.startsWith('#') }.toSet()
                def union_genes = top_genes + bottom_genes
                ali_tuples_from_dir(ali_dir, union_genes)
            }
        }

        // ── Convert + filter: batch or single-gene ───────────────────────────
        def prepBatchSize = (params.selection_prep_batch_size ?: 1) as int
        def filtered_fasta_ch

        if (prepBatchSize > 1) {
            def batches_ch = ali_ch
                .collate(prepBatchSize)
                .map { batch ->
                    def batchID = "prep_batch_${java.util.UUID.randomUUID().toString().replace('-', '').take(12)}"
                    def manifestText = createBatchManifestText(
                        batch.collect { row -> "${row[0]}\t${row[1].name}" }
                    )
                    tuple(batchID, batch.size(), manifestText,
                          batch.collect { row -> row[1] })
                }

            def batched_out = PREP_ALIGNMENTS_BATCHED(batches_ch, tree_val).filtered_fastas
            filtered_fasta_ch = batched_out
                .transpose()
                .map { batchID, f -> tuple(f.name.replace('.filtered.fa', ''), f) }
        } else {
            filtered_fasta_ch = ali_ch
                .map { gid, f -> tuple(gid, f) }
                .combine(tree_val)
                .map { gid, f, tree -> tuple(gid, f, tree) }
            filtered_fasta_ch = PREP_ALIGNMENTS(filtered_fasta_ch).filtered_fasta
        }

        // ── Route filtered FASTAs into top / bottom direction channels ────────
        // multiMap forks the single channel into two independent queues.
        // In gene_set mode, each queue is then filtered to the respective gene set.
        def fanout = filtered_fasta_ch.multiMap { gid, fa ->
            top:    tuple(gid, fa)
            bottom: tuple(gid, fa)
        }

        prep_filtered_top_ch = Channel.empty()
        prep_filtered_bottom_ch = Channel.empty()

        if (mode == 'gene_set') {
            // Keep only genes that appear in the respective directional gene set.
            prep_filtered_top_ch = fanout.top
                .combine(gs_top_val)
                .filter { gid, fa, gsf ->
                    gsf.readLines()
                        .collect { it.trim() }.findAll { it && !it.startsWith('#') }
                        .contains(gid)
                }
                .map { gid, fa, gsf -> tuple(gid, fa) }

            prep_filtered_bottom_ch = fanout.bottom
                .combine(gs_bottom_val)
                .filter { gid, fa, gsf ->
                    gsf.readLines()
                        .collect { it.trim() }.findAll { it && !it.startsWith('#') }
                        .contains(gid)
                }
                .map { gid, fa, gsf -> tuple(gid, fa) }
        } else {
            // 'all' mode: every gene goes to both directions.
            prep_filtered_top_ch    = fanout.top
            prep_filtered_bottom_ch = fanout.bottom
        }

    emit:
        filtered_fasta_top_ch    = prep_filtered_top_ch
        filtered_fasta_bottom_ch = prep_filtered_bottom_ch
        top_species              = top_species_val
        bottom_species           = bottom_species_val
        tree_ch                  = tree_val
}
