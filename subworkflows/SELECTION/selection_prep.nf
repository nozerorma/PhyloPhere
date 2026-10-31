#!/usr/bin/env nextflow

/*
 * SELECTION / shared alignment preprocessing
 *
 * SELECTION_PREP runs the alignment preparation pipeline ONCE and shares the
 * results with FADE:
 *
 *   1. EXTRACT_EXTREME_SPECIES  – top/bottom species lists from the pre-Dunn
 *                                  candidate pool (3.CI-composition.Rmd)
 *   2. [Alignment channel]      – build gene list (all / toy_mode)
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

include { EXTRACT_EXTREME_SPECIES } from './selection_utils.nf'
include { listAlignmentFiles; sampleAlignmentFiles } from '../CT/ct_alignment_files'
include { selectionQuote; selectionHash; selectionSingle; selectionPositive } from './selection_contract'


// ─── Helpers ─────────────────────────────────────────────────────────────────

/**
 * Render a list of TSV row strings into a heredoc-safe manifest block.
 */
def createBatchManifestText(List<String> rows) {
    return rows
        .collect { row -> row.replaceFirst(/^\s+/, '') }
        .join(System.lineSeparator()) + System.lineSeparator()
}

/**
 * Gene ids of the toy sample: n alignments chosen by the seeded shuffle of the name-sorted list that CT draws from
 * (sampleAlignmentFiles), so the same seed picks the same genes whatever order the filesystem lists them in.
 * Returned in the order drawn.
 */
def fade_toy_genes(String ali_dir, n, seed) {
    def files = ali_tuples_from_dir(ali_dir, null).collect { t -> t[1] }
    return sampleAlignmentFiles(files, n, seed).collect { f -> f.baseName }
}

/**
 * Build a List of [gene_id, alignment_File] tuples synchronously (no direction).
 * Returns a plain Groovy List (not a Channel) so it is safe to call inside a
 * channel operator's closure (e.g. flatMap).
 *
 *   wanted  — Set<String> of gene IDs to keep, or null to keep all (empty keeps none).
 */
def ali_tuples_from_dir(String ali_dir, Set wanted) {
    def ali_path = file(ali_dir)
    if (!ali_path.isDirectory()) {
        error "SELECTION_PREP: alignment directory not found: ${ali_dir}"
    }
    // The shared listing is in file-name order and leaves out tables and logs; the extensions below are the ones
    // this module converts (none of them is among those the shared listing leaves out, so they pick the same files).
    def all_files = listAlignmentFiles(ali_dir).findAll { f ->
        def name = f.name.toLowerCase()
        name.endsWith('.phy') ||
        name.endsWith('.phylip') ||
        name.endsWith('.aln') ||
        name.endsWith('.fa') ||
        name.endsWith('.fasta') ||
        !f.name.contains('.')
    }

    if (!all_files) error "SELECTION_PREP: no supported alignment files in ${ali_dir}"
    def ids = all_files.collect { it.baseName }
    if (ids.unique(false).size() != ids.size())
        error "SELECTION_PREP: alignment files must have unique gene IDs (filename without extension)"
    all_files.each { f ->
        def gid = f.baseName
        assert gid && !(gid in ['.', '..']) && gid == gid.trim() &&
            !(f.name =~ /[\t\r\n]/).find() : "SELECTION_PREP: invalid alignment filename: ${f.name}"
    }
    if (wanted != null && !ids.containsAll(wanted)) {
        def message = "SELECTION_PREP: discovery genes have no supported alignment: ${(wanted - ids.toSet()).sort()}"
        log.error message
        error message
    }
    def tuples = all_files.findAll { wanted == null || wanted.contains(it.baseName) }
        .collect { f -> [f.baseName, f] }
    log.info "SELECTION_PREP: ${tuples.size()} genes queued for preprocessing"
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
    cache 'deep'
    publishDir "${params.outdir}/selection/prep_audit", mode: 'copy', overwrite: true, pattern: '*.prep.json'

    input:
    tuple val(gene_id), path(alignment_file), val(alignment_hash)
    tuple path(tree), val(tree_hash)
    val code_hash

    output:
    tuple val(gene_id), path('*.filtered.fa'), emit: filtered_fasta, optional: true
    path '*.prep.json', emit: audit

    script:
    def pyBin = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python' : 'python'
    """
    echo '${code_hash}:${alignment_hash}:${tree_hash}' > inputs.sha256
    ${pyBin} -B ${selectionQuote("${baseDir}/subworkflows/SELECTION/local/src/prepare_alignments.py")} \
        --tree ${selectionQuote(tree)} --gene ${selectionQuote(gene_id)} --alignment ${selectionQuote(alignment_file)}
    """
}

// The same runner handles individual genes and batches. Every input gene gets an
// audit record; technical failures abort, biological exclusions are explicit.
process PREP_ALIGNMENTS_BATCHED {
    tag "${batchID}"
    cache 'deep'
    publishDir "${params.outdir}/selection/prep_audit", mode: 'copy', overwrite: true, pattern: '*.prep.json'

    input:
    tuple val(batchID), val(manifestText), path(alignments, stageAs: 'alignments/*'), val(alignment_hash)
    tuple path(tree), val(tree_hash)
    val code_hash

    output:
    tuple val(batchID), path("*.filtered.fa"), emit: filtered_fastas, optional: true
    path '*.prep.json', emit: audit

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local/src"
    def runnerMode = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'
    """
    echo '${code_hash}:${alignment_hash}:${tree_hash}' > inputs.sha256
    cat <<'MANIFEST_EOF' > ${batchID}.manifest.tsv
${manifestText}MANIFEST_EOF
    bash ${selectionQuote("${local_dir}/run_prep_alignments_batch.sh")} \
        --batch-id ${selectionQuote(batchID)} --manifest ${batchID}.manifest.tsv \
        --workers ${task.cpus} --runner-mode ${runnerMode} \
        --tree ${selectionQuote(tree)} --local-dir ${selectionQuote(local_dir)}
    """
}


// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow SELECTION_PREP {

    take:
        species_source_input // Channel<path> to the pre-Dunn candidate fg/bg pool
                             // (candidate_species.tab, species\tcontrast_group\tpair) from
                             // 3.CI-composition.Rmd; or null → falls back to canonical outdir path
        tree_input         // Channel<path> or null → falls back to params.tree
        pp_top_ch          // Legacy input, unused (kept for backwards-compatible signature)
        pp_bottom_ch       // Legacy input, unused (kept for backwards-compatible signature)
        ct_discovery_input // Optional CT discovery table for toy_mode gene reuse

    main:

        def sourceNames = ['selection_prep.nf', 'selection_utils.nf', 'selection_contract.nf',
            'local/src/alignment_io.py', 'local/src/prepare_alignments.py',
            'local/src/extract_extreme_species.py', 'local/src/run_prep_alignments_batch.sh']
        def code_hash = selectionHash(sourceNames.collect { file("${baseDir}/subworkflows/SELECTION/${it}") })
        def prepBatchSize = selectionPositive(params.selection_prep_batch_size != null ? params.selection_prep_batch_size : 1, 'selection_prep_batch_size')
        def toyN = params.toy_mode ? selectionPositive(params.toy_n != null ? params.toy_n : 50, 'toy_n') : null
        def toySeed = params.seed != null ? params.seed : 1998
        if (params.toy_mode) assert toySeed.toString() ==~ /-?[0-9]+/ : 'SELECTION_PREP: seed must be an integer'

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
            .map { files -> selectionSingle(files, 'tree') }

        def species_source_ch = (species_source_input ?: Channel.empty())
            .ifEmpty {
                // Pre-Dunn candidate fg/bg pool from CI_COMPOSITION_REPORT (3.CI-composition.Rmd).
                def fallback = "${params.outdir}/data_exploration/1.Data-exploration/5.CI_overlaps/candidate_species.tab"
                def f = file(fallback)
                // Don't assert existence here — file is created by upstream CI_COMPOSITION_REPORT process
                f
            }
            .collect()
            .map { files -> selectionSingle(files, 'species pool') }

        def ali_dir = params.alignment ?:
            error("SELECTION_PREP: no alignment directory specified (--alignment)")

        // ── Extract foreground species (once, shared with FADE) ───
        def species_lists      = EXTRACT_EXTREME_SPECIES(species_source_ch.map { f -> tuple(f, selectionHash([f])) }, code_hash)
        top_species_val        = species_lists.top_species.collect().map { it[0] }
        bottom_species_val     = species_lists.bottom_species.collect().map { it[0] }

        // ── Build alignment channel (gene_id, alignment_file) ────────────────
        // Conversion and tree-filtering are direction-agnostic;
        // direction is introduced later inside FADE.
        def ali_ch

        if (params.toy_mode) {
            def resolved_ct = (ct_discovery_input ?: Channel.empty())
                .ifEmpty { file('NO_FILE') }
                .collect()
                .map { inputs ->
                    assert inputs.size() == 1 : 'SELECTION_PREP: expected at most one CT discovery table'
                    inputs[0]
                }

            ali_ch = resolved_ct.flatMap { discovery_file ->
                def toy_genes = null
                if (discovery_file.name != 'NO_FILE') {
                    assert discovery_file.isFile() : "SELECTION_PREP: discovery not found: ${discovery_file}"
                    def lines = discovery_file.readLines()
                    def fields = lines ? lines[0].split('\t', -1).toList() : []
                    assert fields.count('gene') == 1 : 'SELECTION_PREP: discovery needs a unique gene column'
                    def geneColumn = fields.indexOf('gene')
                    toy_genes = lines.drop(1)
                        .collect { line ->
                            def row = line.split('\t', -1)
                            assert row.size() == fields.size() && row[geneColumn] : 'SELECTION_PREP: malformed discovery row'
                            row[geneColumn]
                        }
                        .findAll { it }
                        .toSet()
                    log.info "[toy_mode] SELECTION_PREP: reusing ${toy_genes.size()} genes from CT discovery"
                } else {
                    // Seeded shuffle of the name-sorted list: reproducible, and the same genes CT draws from this directory.
                    toy_genes = fade_toy_genes(ali_dir, toyN, toySeed as long).toSet()
                    log.info "[toy_mode] SELECTION_PREP: using ${toy_genes.size()} randomly sampled genes (seed=${toySeed})"
                }
                ali_tuples_from_dir(ali_dir, toy_genes)
            }
        } else {
            ali_ch = Channel.fromList(ali_tuples_from_dir(ali_dir, null))
        }

        // ── Convert + filter: batch or single-gene ───────────────────────────
        def tree_fingerprinted = tree_val.map { f -> tuple(f, selectionHash([f])) }
        def filtered_fasta_ch
        def prep_audit_ch

        if (prepBatchSize > 1) {
            def prepBatchCounter = 0
            def batches_ch = ali_ch
                .toSortedList({ a, b -> a[0] <=> b[0] })
                .flatMap()
                .collate(prepBatchSize)
                .map { batch ->
                    def batchID = String.format('prep_batch_%05d', ++prepBatchCounter)
                    def manifestText = createBatchManifestText(
                        batch.collect { row -> "${row[0]}\t${row[1].name}" }
                    )
                    def paths = batch.collect { row -> row[1] }
                    tuple(batchID, manifestText, paths, selectionHash(paths))
                }

            def batched_out = PREP_ALIGNMENTS_BATCHED(batches_ch, tree_fingerprinted, code_hash)
            prep_audit_ch = batched_out.audit.flatten()
            filtered_fasta_ch = batched_out.filtered_fastas
                .transpose()
                .map { batchID, f -> tuple(f.name.replaceFirst(/\.filtered\.fa$/, ''), f) }
        } else {
            filtered_fasta_ch = ali_ch
                .map { gid, f -> tuple(gid, f, selectionHash([f])) }
            def single_out = PREP_ALIGNMENTS(filtered_fasta_ch, tree_fingerprinted, code_hash)
            filtered_fasta_ch = single_out.filtered_fasta
            prep_audit_ch = single_out.audit
        }

        // ── Route filtered FASTAs into top / bottom direction channels ────────
        // multiMap forks the preprocessed gene fastas into both directional queues.
        def fanout = filtered_fasta_ch.multiMap { gid, fa ->
            top:    tuple(gid, fa)
            bottom: tuple(gid, fa)
        }

        def prep_filtered_top_ch    = fanout.top
        def prep_filtered_bottom_ch = fanout.bottom

    emit:
        filtered_fasta_top_ch    = prep_filtered_top_ch
        filtered_fasta_bottom_ch = prep_filtered_bottom_ch
        top_species              = top_species_val
        bottom_species           = bottom_species_val
        tree_ch                  = tree_val
        prep_audit               = prep_audit_ch
}
