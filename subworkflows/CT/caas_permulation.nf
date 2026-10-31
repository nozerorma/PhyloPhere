#!/usr/bin/env nextflow

// caas_permulation.nf — Permulation core: replay of the permulated labelings and the observed one, and the genome-wide null.
// PhyloPhere | subworkflows/CT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CAAS permulation-excess null: builds a genome-wide excess null for the CAAS FCS
 *  pathway enrichment, and scores the observed labeling (b_0) with the same code.
 *
 *    1. SUBSET_RESAMPLE_PERMS takes the first N (caas_full_perms) permulated labelings
 *       in cycle order (N = 0 keeps none) and writes them in blocks of caas_cycle_block_size
 *       consecutive cycles (blocks/blk_NNNN). The real labeling b_0 is rebuilt from the
 *       observed design into its own directory (blocks/b0); it is never part of the null
 *       (its shards go to caas_permulation/b0/).
 *    2. CAAS_CORE_BATCHED (one task per tile: a batch of ct_core_batch_size genes against
 *       one block, or against b_0) replays the block's labelings over each gene's alignment and
 *       its ASR in process (disambiguation_perms_main.py) and writes the tile tables of the
 *       null: per (gene, cycle) a unit row and the scored positions of the cycle (null_tables/,
 *       Parquet). The b_0 tile also replays b_0 through `ct perm-replay` and scores it as full
 *       records (observed_b0_main.py) into b0_observed/.
 *    3. CAAS_CORE_OBSERVED writes the observed contract files (discovery.tab,
 *       background.output, background_genes.output, meta_caas/,
 *       ct_disambiguation/caas_convergence_master.csv) from the b_0 tiles' slices.
 *    4. CAAS_CORE_MERGE merges the tiles of one block of cycles (all the genes of a cycle are in
 *       the same block, which gene removal and the reference pools of the gene score need):
 *       gene removal, the pools of the cycle, the gene statistics as genes x cycles matrices,
 *       and the scored positions of the block. CAAS_CORE_MERGE_REDUCE joins the blocks into
 *       the validated block index, roster and small tables. A block merge needs the gene universe, which
 *       the post-processing of the observed results provides, so it runs after CAAS_CORE_OBSERVED.
 *
 *  The canonical blocks feed the FCS p.perm path, which gives the CAAS
 *  scoring FCS report a permulation-corrected p.perm, as RER does.
 *
 *  Consumes:  resample directory (RESAMPLE), observed design (trait file or
 *             traitfile_H*.tab directory), alignments, species tree
 *  Produces:  caas_permulation/ (resample_perms.tab, optional caas_perms.rds, cycle_roster.txt,
 *             blocks/<block>/, b0/, ...) and the observed contract files
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


include { observedSourceHash as scoringSourceHash; observedQuote as scoringQuote; observedInputHash as scoringInputHash } from '../CT_DISAMBIGUATION/ct_observed'
def core_scoring_source_hash = scoringSourceHash("${baseDir}/subworkflows/CT_DISAMBIGUATION/local", "${baseDir}/subworkflows/CT/local/modules/design_io.py", "${baseDir}/subworkflows/CT/local") + java.security.MessageDigest.getInstance('SHA-256').digest(file("${baseDir}/subworkflows/SCORING/local/src/assemble_block_matrices.R").bytes).encodeHex().toString()

// ── Helpers ──────────────────────────────────────────────────────────────────

// True unless params.caas_perms_postproc is off (Boolean false, or the strings false / 0 / no).
def caasPostprocOn() {
    def raw = params.containsKey('caas_perms_postproc') ? params.caas_perms_postproc : true
    return (raw instanceof Boolean) ? raw : !(raw?.toString()?.toLowerCase() in ['false', '0', 'no'])
}

// True when the core also writes the strict tables: params.strict_conservation_output, and a discovery that tolerates
// conserved pairs (min_divergent_fraction < 1); with fraction 1 the strict rows are all the rows.
def strictOutputOn() {
    return (params.strict_conservation_output as boolean) && ((params.min_divergent_fraction ?: 0.5) as double) < 1.0
}

// CLI arguments that give disambiguation_perms_main.py the observed CT_POSTPROC cluster filter for the per-cycle null:
// the trains of every cycle are flagged (they count for the dubious-gene test) and, unless params.remove_caas_clusters is
// off, left out of the scores. Empty when the filter is off or the gene annotation file is a NO_* sentinel. With
// params.caas_map_dir the trains measure their span in untrimmed columns, as the observed CT_FILTER does.
def caasTileArgs(gene_lengths) {
    if (!caasPostprocOn() || new File(gene_lengths.toString()).name =~ /^NO_/) return ""
    def map_arg = params.caas_map_dir ? "--train-map-dir ${scoringQuote(params.caas_map_dir)}" : ""
    def count_arg = params.filter_min_positions != null ? "--clust-min-positions ${scoringQuote(params.filter_min_positions)}" : ""
    return "--postproc-filter ${count_arg} --clust-minlen ${scoringQuote(params.filter_minlen)} --clust-maxcaas ${scoringQuote(params.filter_maxcaas)} ${params.remove_caas_clusters ? '' : '--keep-clusters'} ${map_arg}"
}

// CLI arguments that give merge_null_block.py the observed gene removal (dubious and extreme genes, or the length-stratified
// fences, each cycle against its own units). Empty when the filters are off or the gene annotation file is a NO_* sentinel.
def caasRemovalArgs(gene_lengths) {
    if (!caasPostprocOn() || new File(gene_lengths.toString()).name =~ /^NO_/) return ""
    def map_arg = params.caas_map_dir ? "--map-dir ${scoringQuote(params.caas_map_dir)}" : ""
    return "--gene-lengths ${scoringQuote(gene_lengths)} ${map_arg} --gene-filter-mode ${scoringQuote(params.gene_filter_mode)} --gene-strata ${scoringQuote(params.gene_filter_strata)} --gene-k-far ${scoringQuote(params.gene_filter_k_far)} --gene-k-near ${scoringQuote(params.gene_filter_k_near)}"
}

// Shell lines that write .ct_perm_replay_batch_args, the effective CAAS limits of `ct perm-replay` (patterns, miss_pair,
// max_conserved from min_divergent_fraction and the pair count, min_non_overlapping, gap and missing limits). The replay
// of the permulations and the positional universe read the same file, so the universe is judged with the limits of CT.
def replayArgsScript(caas_config) {
    return """
    # With multi_hypothesis the design is a directory of traitfile_H*.tab (same K): one .tab gives the pair count
    if [ -d ${scoringQuote(caas_config)} ]; then
        _cfg_file=\$(find -L ${scoringQuote(caas_config)} -type f -name '*.tab' | sort | head -n 1)
    else
        _cfg_file=${scoringQuote(caas_config)}
    fi
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$_cfg_file" | sort -nu | wc -l | tr -d ' ')
    _frac() { awk -v n="\$n_pairs" -v f="\$1" 'BEGIN{printf "%d", int(n*f)}'; }
    declare -a extra_opts=(--patterns "${params.patterns}")
    if [ "${params.miss_pair}" = "true" ]; then extra_opts+=(--miss_pair); fi
    extra_opts+=(--max_conserved \$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}'))
    extra_opts+=(--min_non_overlapping "${params.min_non_overlapping ?: 2}")
    extra_opts+=(--max_bg_gaps \$(_frac ${params.max_bg_gaps_fraction}) --max_fg_gaps \$(_frac ${params.max_fg_gaps_fraction}) --max_gaps \$(_frac ${params.max_gaps_fraction}))
    extra_opts+=(--max_bg_miss \$(_frac ${params.max_bg_miss_fraction}) --max_fg_miss \$(_frac ${params.max_fg_miss_fraction}) --max_miss \$(_frac ${params.max_miss_fraction}))
    echo "\${extra_opts[@]}" > .ct_perm_replay_batch_args
"""
}

// ── 1. First N permulated labelings in cycle order, plus b_0, split into cycle blocks ─────────
// resample_perms.tab holds b_0 and the N cycles (the null's cycle roster). The cycles are also written in blocks of
// block_size consecutive cycles (the first block_size cycles in cycle order, then the next, ...), so a replay task reads
// only the labelings of its block:
//   blocks/blk_NNNN/resample_perms.tab  the labelings of the block, in cycle order
//   blocks/blk_NNNN/fop_pairs.tsv       the PSS weights of the block's cycles (FOP mirror only)
//   blocks/b0/                          the same two files for the real labeling b_0, never part of a block
process SUBSET_RESAMPLE_PERMS {
    tag "caas_perms_subset|N=${n_perms}|block=${block_size}"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: '{resample_perms.tab,ct_preparation.json}'

    input:
    path resample_dir
    val  n_perms
    val  block_size
    path caas_config   // observed design: trait file or multi-hypothesis dir (source of b_0)

    output:
    path "resample_perms.tab", emit: subset
    path "ct_preparation.json", emit: provenance
    path "blocks/*", emit: blocks   // blocks/b0 and one blocks/blk_NNNN per block of cycles

    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    ${run} "${baseDir}/subworkflows/CT/local/scripts/prepare_labelings.py" \
        --source "${resample_dir}" --config "${caas_config}" \
        --n "${n_perms}" --block-size "${block_size}" \
        --multi "${params.multi_hypothesis == true}" \
        --allow-relaxed "${params.ct_allow_relaxed_null == true}"
    """
}

// ── 2a. Positional universe of the observed design, one task per gene batch ─────────────────────
// The positions a gene can score under the observed hypotheses H1..Hn: those b_0 tested that pass the CAAS test, with the
// effective limits of CT, under some orientation of the pairs of some hypothesis (`ct perm-replay --export_b0_universe`,
// modules/universe.py; the same kernel as the replay, no pair is added). One `<gene>.b0.universe` line per gene, NULL when
// it has none. Every tile of the batch (the permuted blocks and b_0) reads it, so the null and the observed labeling are
// scored over the same positions, and CAAS_CORE_OBSERVED writes it as background.output. With the strict tables on, the
// strict universe (`<gene>.b0.universe_strict`) is the same test with no conserved pair allowed.
process CAAS_UNIVERSE {
    tag "$batchID (${batchSize} genes)"
    label 'process_resample'

    input:
    tuple val(batchID), val(batchSize), val(manifestText), path(alignmentFiles, stageAs: 'alignments/*')
    path b0_block, stageAs: 'b0_block'   // blocks/b0 of SUBSET_RESAMPLE_PERMS: the b_0 hypotheses
    file caas_config

    output:
    tuple val(batchID), path("universe"), emit: universe

    script:
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh ${scoringQuote("${baseDir}/subworkflows/CT/local/ct")}"
        : scoringQuote("${baseDir}/subworkflows/CT/local/ct")
    """
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cat > ${batchID}.manifest.tsv <<'EOF'
""" + manifestText + """EOF
    ${replayArgsScript(caas_config)}
    bash ${scoringQuote("${baseDir}/subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh")} \\
        --batch-id ${batchID} \\
        --manifest ${batchID}.manifest.tsv \\
        --caas-config ${scoringQuote(caas_config)} \\
        --resampled-path ${scoringQuote("b0_block/resample_perms.tab")} \\
        --workers ${task.cpus} \\
        --ali-format ${scoringQuote(params.ali_format)} \\
        --ct-bin ${scoringQuote(ctBinary)} \\
        --export-universe 1 ${strictOutputOn() ? '--export-universe-strict 1' : ''} \\
        --extra-args-file .ct_perm_replay_batch_args
    mkdir universe
    find . -maxdepth 1 \\( -name '*.b0.universe' -o -name '*.b0.universe_strict' \\) -exec mv {} universe/ \\;
    """
}

// ── 2. One task per (gene batch, cycle block): the replay of the block's labelings over the batch's genes ───
// A tile is the replay of one batch of genes against the labelings of one cycle block (or of b_0, the tile
// whose block directory is blocks/b0). disambiguation_perms_main.py replays the block's labelings over
// each gene's alignment in process (the CAAS kernel of `ct perm-replay`, one base cycle at a time), scores each
// cycle against the gene's ASR and reduces it to a (gene, cycle) unit row and the scored positions of the cycle; no
// discovery export is written, and a gene is replayed in one worker, so its ASR context is loaded once per tile.
// The tables of the batch's genes are one pair of Parquet files (src/utils/null_tables.py). What needs every gene of a
// cycle (gene removal, the pools of the gene score) is done once per block, in CAAS_CORE_MERGE.
// The b_0 tile also runs `ct perm-replay` for the b_0 slice of each gene (the discovery.tab rows of the real labeling
// and the positions it tested); observed_b0_main.py then scores those rows as full records.
// Output: null_tables/ (the tile's tables, with the name of its block) and, for the b_0 tile, b0_observed/ (per gene
// with b_0 hits: its discovery rows, its tested positions and its master rows).
process CAAS_CORE_BATCHED {
    ext scoring_source_hash: core_scoring_source_hash
    tag "$batchID (${batchSize} genes) ${block.name}"
    label 'process_resample'

    input:
    tuple val(batchID), val(batchSize), val(manifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(block), path(universe, stageAs: 'universe_in')   // block: blocks/blk_NNNN or blocks/b0 of SUBSET_RESAMPLE_PERMS; universe: the batch's CAAS_UNIVERSE output
    file caas_config
    path tree_file
    path gene_lengths, stageAs: 'inputs/gene_lengths/*' // gene_ensembl_file (CT_POSTPROC filters) or NO_FILE
    path taxid_map, stageAs: 'inputs/taxid/*' // tax_id map of the species (curated by NAME_CURATION, or params.tax_id) or the NO_FILE sentinel

    output:
    tuple val("${block.name}"), path("null_tables"), emit: null_tables                         // the tile tables of the block
    tuple val("${block.name}"), path("null_tables_strict"), emit: null_tables_strict, optional: true   // strict pooling of the same replay (strictOutputOn)
    path "b0_observed", emit: b0_observed, optional: true                    // b_0 tile only
    path "b0_observed_strict", emit: b0_observed_strict, optional: true      // b_0 rows without conserved pairs, scored on their own (b_0 tile only)

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh ${scoringQuote("${baseDir}/subworkflows/CT/local/ct")}"
        : scoringQuote("${baseDir}/subworkflows/CT/local/ct")
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = new File(taxid_map.toString()).name != 'NO_FILE' ? taxid_map : ''
    def ensembl_file = new File(gene_lengths.toString()).name.startsWith('NO_') ? '' : gene_lengths
    def max_tasks_per_child = params.ct_disambig_max_tasks_per_child ?: 50
    def postproc_maps_hash = params.caas_map_dir ? scoringInputHash(params.caas_map_dir) : ''
    def tile_args = caasTileArgs(gene_lengths)
    // The strict tables come from the same replay: the null tables pooled without the hypotheses that have conserved
    // pairs, and the b_0 rows without conserved pairs scored again (the observed labeling is a single labeling).
    def strict_on = strictOutputOn()
    // The b_0 tile replays the real labeling alone; the block tiles hold permuted cycles only.
    def is_b0 = (block.name == 'b0')
    def strict_flag = strict_on ? '--strict-output' : ''
    def strict_observed_flag = strict_on ? '--strict-output-dir b0_observed_strict' : ''
    // The observed labeling as full records: one master shard per gene with b_0 hits (the b_0 tile only).
    def observed_b0 = !is_b0 ? '' : """
    # The observed labeling as full records: one master shard per gene with b_0 hits.
        ${run} ./observed_b0_main.py \\
            --alignment-dir alignments \\
            --tree ${scoringQuote(tree_file)} \\
            --b0-dir b0_observed ${strict_observed_flag} \\
            --alignment-manifest ${batchID}.manifest.tsv --manifest-alignment-dir alignments \\
            --ali-format ${scoringQuote(params.ali_format)} \\
            --design ${scoringQuote(caas_config)} \\
            \$FOP_PAIRS_ARG \\
            --output-dir b0_observed \\
            --asr-model ${scoringQuote(params.ct_disambig_asr_model)} \\
            --posterior-threshold ${scoringQuote(params.ct_disambig_posterior_threshold)} \\
            --workers ${task.cpus} \\
            --max-tasks-per-child ${max_tasks_per_child} \\
            --asr-cache-dir ${scoringQuote(asr_cache_dir)} \\
            ${taxid_mapping ? "--taxid-mapping ${scoringQuote(taxid_mapping)}" : ''} \\
            ${ensembl_file ? "--ensembl-genes-file ${scoringQuote(ensembl_file)}" : ''}
    """
    // The b_0 slice of each gene through `ct perm-replay` (the b_0 tile only): discovery.tab rows and tested positions.
    def ct_b0 = !is_b0 ? '' : """bash ${scoringQuote("${baseDir}/subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh")} \\
        --batch-id ${batchID} \\
        --manifest ${batchID}.manifest.tsv \\
        --caas-config ${scoringQuote(caas_config)} \\
        --resampled-path ${scoringQuote("${block}/resample_perms.tab")} \\
        --workers ${task.cpus} \\
        --ali-format ${scoringQuote(params.ali_format)} \\
        --ct-bin ${scoringQuote(ctBinary)} \\
        --export-groups 0 \\
        --export-perm-discovery 0 \\
        --export-b0 1 \\
        --extra-args-file .ct_perm_replay_batch_args
    find . -maxdepth 1 \\( -name '*.b0.discovery.tsv' -o -name '*.b0.background' \\) -exec mv {} b0_observed/ \\;
    # the positional universe of the batch's genes goes with the b_0 slice: background.output is written from it
    cp universe_in/*.b0.universe* b0_observed/"""
    def strict_move = strict_on ? "mv caas_perms_out/null_tables_strict null_tables_strict" : ''
    """
    # CT scoring contract v4: validated provenance, strict overlap, complete observed/null scoring.
    echo '${task.ext.scoring_source_hash}:${postproc_maps_hash}' > scoring_source.sha256
    # One worker pool runs the genes in parallel; BLAS and OpenMP stay at one thread each.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    ${is_b0 ? 'mkdir -p b0_observed' : ''}
    ${is_b0 && strict_on ? 'mkdir -p b0_observed_strict' : ''}
    cat > ${batchID}.manifest.tsv <<'EOF'
""" + manifestText + """EOF

    ${replayArgsScript(caas_config)}

    # The b_0 slice of each gene: the discovery.tab rows of the real labeling and the positions it tested.
    ${ct_b0}

    # The PSS weights of the block's cycles (absent for a design without FOP hypotheses).
    FOP_PAIRS_ARG=""
    if [ -f ${block}/fop_pairs.tsv ]; then FOP_PAIRS_ARG="--fop-pairs ${block}/fop_pairs.tsv"; fi

    cp -R ${scoringQuote(local_dir)}/. .
    cp ${scoringQuote("${baseDir}/subworkflows/CT/local/modules/design_io.py")} src/core/_ct_design_io.py
    find . -name '*.pyc' -delete 2>/dev/null || true
    mkdir -p caas_perms_out
    ${run} ./disambiguation_perms_main.py \\
        --alignment-dir alignments \\
        --tree ${scoringQuote(tree_file)} \\
        --manifest ${batchID}.manifest.tsv --manifest-alignment-dir alignments \\
        --design ${scoringQuote(caas_config)} \\
        --ct-dir ${scoringQuote("${baseDir}/subworkflows/CT/local")} \\
        --replay-args-file .ct_perm_replay_batch_args \\
        --universe-dir universe_in \\
        --ali-format ${scoringQuote(params.ali_format)} \\
        --resample-dir ${block} \\
        --output-dir caas_perms_out \\
        ${strict_flag} \\
        \$FOP_PAIRS_ARG ${tile_args} \\
        --asr-model ${scoringQuote(params.ct_disambig_asr_model)} \\
        --posterior-threshold ${scoringQuote(params.ct_disambig_posterior_threshold)} \\
        --workers ${task.cpus} \\
        --max-tasks-per-child ${max_tasks_per_child} \\
        --asr-cache-dir ${scoringQuote(asr_cache_dir)} \\
        ${taxid_mapping ? "--taxid-mapping ${scoringQuote(taxid_mapping)}" : ''} \\
        ${ensembl_file ? "--ensembl-genes-file ${scoringQuote(ensembl_file)}" : ''}
    mv caas_perms_out/null_tables null_tables
    ${strict_move}
    rm -rf caas_perms_out

    ${observed_b0}
    """
}

// ── Permulation core workflow ────────────────────────────────────────────────
// Tiles: batches of `ct_core_batch_size` genes (1 = one task per gene), in gene-name order, crossed with the cycle
// blocks of SUBSET_RESAMPLE_PERMS and its b_0 directory.
workflow CAAS_CORE {
    take:
        align_tuple      // Channel<tuple(id, alignmentFile)> to replay
        caas_config
        resample_subset  // resample_perms.tab: b_0 and the N cycles (the null's cycle roster)
        blocks           // Channel of the block directories of SUBSET_RESAMPLE_PERMS (blk_NNNN and b0), as one list or one by one
        tree_file        // species tree
        gene_lengths
        taxid_map        // tax_id map of the species (curated by NAME_CURATION, or params.tax_id) or NO_FILE

    main:
        def batchSize = (params.ct_core_batch_size ?: 20) as int
        def liveCounter = 0
        def live = align_tuple
            .toSortedList({ a, b -> a[0] <=> b[0] })
            .flatMap()
            .collate(batchSize)
            .map { batch ->
                def id = String.format('caas_core_batch_%05d', ++liveCounter)
                def manifest = batch.collect { row -> "${row[0]}\t${row[1].name}" }.join('\n') + '\n'
                tuple(id, batch.size(), manifest, batch.collect { row -> row[1] }.unique { f -> f.name })
            }
        // The positional universe of each gene batch (from b_0 and the observed pairs) is computed once and read by every
        // tile of the batch, the permuted blocks and b_0 alike.
        def b0_block = blocks.flatten().filter { it.name == 'b0' }.collect().map { items -> items[0] }
        def universe_cfg = caas_config.collect().map { items -> items[0] }
        def universe_out = CAAS_UNIVERSE(live, b0_block, universe_cfg)
        // Every gene batch replays against every block directory.
        def tiles = live.combine(blocks.flatten()).combine(universe_out.universe, by: 0)
        // The single-item inputs (resample subset, gene lengths, trait file, tree) are queue
        // channels once they have crossed a take:/emit: boundary, not value channels. A process pairs its
        // input channels positionally and stops when any one runs out, so paired with the many-item tile
        // channel each of them would end CAAS_CORE_BATCHED after its first tile. .collect().map { items -> items[0] }
        // turns each into a value channel that every tile can read.
        def subset_bc = resample_subset.collect().map { items -> items[0] }
        def lengths_bc = gene_lengths.collect().map { items -> items[0] }
        def config_bc = caas_config.collect().map { items -> items[0] }
        def tree_bc   = tree_file.collect().map { items -> items[0] }
        def taxid_bc  = taxid_map.collect().map { items -> items[0] }
        def core = CAAS_CORE_BATCHED(tiles, config_bc, tree_bc, lengths_bc, taxid_bc)

    emit:
        labelings      = subset_bc                 // the labelings file the tiles replayed (the null's cycle roster)
        null_tables    = core.null_tables          // tuple(block name, tile tables), one per tile
        null_tables_strict = core.null_tables_strict   // empty unless the strict tables are on
        blocks         = blocks.flatten().map { b -> tuple(b.name, b) }   // tuple(block name, block directory)
        universe       = universe_out.universe     // tuple(batch id, directory of its *.b0.universe files)
        b0_observed    = core.b0_observed
        b0_observed_strict = core.b0_observed_strict
}

// ── 3. The observed contract files, from the b_0 slices of the batches ───────
// contract_main.py turns the b0_observed/ directories into discovery.tab, background.output,
// background_genes.output, meta_caas/ and ct_disambiguation/caas_convergence_master.csv: every
// table holds the genes in name order, whatever order the batches finished in. Batches that
// reuse perm-discovery exports carry no b_0 slice and give no file. The files are published
// where the rest of the pipeline reads them: caastools/, meta_caas/meta_caas/ and
// ct_disambiguation/.
process CAAS_CORE_OBSERVED {
    ext scoring_source_hash: core_scoring_source_hash
    tag "caas_core_observed"
    label 'process_low'
    // ext.publish_prefix is empty for this process and "/strict" for its strict alias (CAAS_CORE_OBSERVED_STRICT), which
    // writes the same files from the strict b_0 directories.
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}/caastools" }, mode: 'copy', overwrite: true, pattern: '{discovery.tab,background.output,background_genes.output}'
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}/meta_caas" }, mode: 'copy', overwrite: true, pattern: 'meta_caas'
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}" }, mode: 'copy', overwrite: true, pattern: 'ct_disambiguation'

    input:
    path b0Observed, stageAs: 'b0obs_*'   // the batches' b0_observed directories
    path design                           // observed design: the master columns come from it

    output:
    path "discovery.tab",                emit: discovery, optional: true
    path "background.output",            emit: background, optional: true
    path "background_genes.output",      emit: background_genes, optional: true
    path "ct_disambiguation",            emit: results_dir, optional: true
    path "ct_disambiguation/caas_convergence_master.csv", emit: master_csv, optional: true
    path "meta_caas",                    emit: meta_caas_dir, optional: true
    path "meta_caas/meta_caas.tsv", emit: meta_caas, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    echo '${task.ext.scoring_source_hash}' > scoring_source.sha256
    cp -R ${scoringQuote(local_dir)}/. .
    cp ${scoringQuote("${baseDir}/subworkflows/CT/local/modules/design_io.py")} src/core/_ct_design_io.py
    find . -name '*.pyc' -delete 2>/dev/null || true
    ${py} ./contract_main.py --b0-dirs b0obs_* --design ${scoringQuote(design)} --output-dir .
    """
}

// ── 4. Merge the tiles of each block of cycles and join the blocks ───────────────
// CAAS_CORE_MERGE runs once per block (b_0 is its own block). It merges the
// cached replay tiles, applies gene removal within each cycle and writes raw
// maximum gene matrices plus positive-position counts. No ECDF pool or exponent
// is computed. Observed and null gene scores both use max; FCS calibrates the
// same statistic on both sides. Position-level null products are unchanged.
process CAAS_CORE_MERGE {
    ext scoring_source_hash: core_scoring_source_hash
    cache 'deep'
    tag "caas_core_merge|${block}"
    label 'process_medium'
    // ext.publish_prefix is "/strict" for the strict alias (CAAS_CORE_MERGE_STRICT)
    // the output directory is named after the block: caas_permulation/b0/ and caas_permulation/blocks/<block>/
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}/caas_permulation${block == 'b0' ? '' : '/blocks'}" }, mode: 'copy', overwrite: true

    input:
    tuple val(block), path(tiles, stageAs: 'tile_*'), path(block_dir, stageAs: 'block_input')   // the block's tile tables and its directory of SUBSET_RESAMPLE_PERMS (staged under a name other than the output's)
    path universe       // gene universe (cleaned_background_main.txt) or NO_FILE
    path gene_lengths   // gene_ensembl_file (gene removal) or NO_FILE

    output:
    tuple val(block), path("${block}"), emit: block_out

    script:
    def postproc_maps_hash = params.caas_map_dir ? scoringInputHash(params.caas_map_dir) : ''
    def disambig_local = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scoring_local  = "${baseDir}/subworkflows/SCORING/local"
    // The sentinel convention is a NO_* basename (NO_FILE here, NO_BACKGROUND from SCORING's resolved_background), so
    // the prefix is matched, not one literal.
    def universe_arg   = universe.name.startsWith('NO_') ? "" : "--universe ${universe}"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript'  : 'Rscript'
    def b0_arg = block == 'b0' ? '--include-b0' : ''
    """
    echo '${task.ext.scoring_source_hash}:${postproc_maps_hash}' > scoring_source.sha256
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    # merge_null_block.py imports its src package relative to its own location.
    cp -R ${disambig_local}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    ${py} ./merge_null_block.py \
        --tiles tile_*/ \
        --labelings ${block_dir}/resample_perms.tab \
        --output-dir ${block} ${b0_arg} ${universe_arg} \
        --gene-stat ${params.gene_stat ?: 'max'} \
        --seed ${params.seed ?: 1998} ${caasRemovalArgs(gene_lengths)}

    # The matrices go to one RDS per block; the raw files are only the hand-over.
    ${rs} ${scoring_local}/src/assemble_block_matrices.R --dir ${block}/matrices --output ${block}/matrices.rds
    rm -rf ${block}/matrices
    """
}

// The blocks joined: the roster (N of every consumer of the null), the small tables of the whole null, the positions the
// null scores (nd per position, the universe of the accumulation test and of the position p).
// Canonical matrices stay in blocks; caas_export_single requests a compatibility export.
process CAAS_CORE_MERGE_REDUCE {
    tag "caas_core_merge_reduce"
    label 'process_medium'
    // ext.publish_prefix is "/strict" for the strict alias (CAAS_CORE_MERGE_REDUCE_STRICT)
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}/caas_permulation" }, mode: 'copy', overwrite: true,
               pattern: '{null_metadata.json,caas_perms.rds,perm_pos_sample.tsv,perm_pos_quantiles.tsv,cycle_roster.txt,removed_units.tsv,cycle_summary.tsv,cycle_hist.tsv,null_positions.parquet,counts.json}'

    input:
    path blocks, stageAs: 'blk_*'    // the block directories of CAAS_CORE_MERGE (permuted cycles only)
    path universe                    // gene universe (cleaned_background_main.txt) or NO_FILE
    path labelings                   // complete replay input: detects missing whole blocks

    output:
    path "caas_perms.rds",            emit: perms, optional: true
    path "null_metadata.json",        emit: metadata
    path "removed_units.tsv",         emit: removed_units, optional: true   // null gene removal (only when a gene annotation was given)
    path "cycle_roster.txt",          emit: cycle_roster
    path "perm_pos_sample.tsv",       emit: pos_sample
    path "perm_pos_quantiles.tsv",    emit: pos_quantiles
    path "cycle_summary.tsv",         emit: cycle_summary
    path "cycle_hist.tsv",            emit: cycle_hist
    path "null_positions.parquet",    emit: null_positions
    path "counts.json",               emit: counts

    script:
    def disambig_local = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scoring_local  = "${baseDir}/subworkflows/SCORING/local"
    def universe_arg   = universe.name.startsWith('NO_') ? "" : "--universe ${universe}"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript'  : 'Rscript'
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cp -R ${disambig_local}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    ${py} ./reduce_null_blocks.py --blocks blk_*/ --labelings ${labelings} --output-dir .
    ${rs} ${scoring_local}/src/scoring_caas_perms.R \
        --blocks blk_*/matrices.rds \
        --cycles cycle_roster.txt \
        ${universe_arg} \
        --metadata null_metadata.json \\
        ${params.caas_export_single == true ? '--output caas_perms.rds' : '--validate-only'}
    """
}

// ── Resample subset workflow ─────────────────────────────────────────────────
// Called by workflows/ct.nf, where the resample directory is available, and by main.nf for a run that
// starts from a resample file (--resample_from).
workflow CAAS_PERMS_PREP {
    take:
        caas_config        // observed design (trait file, or directory of traitfile_H*.tab)
        resample_dir       // resample_*.tab directory (RESAMPLE output)

    main:
        def subset = SUBSET_RESAMPLE_PERMS(resample_dir, params.caas_full_perms == null ? 10 : (params.caas_full_perms as int),   // 0 is a valid N: b_0 only
                                           (params.caas_cycle_block_size == null ? 500 : (params.caas_cycle_block_size as int)), caas_config)  // 0 is refused by the process

    emit:
        resample_subset = subset.subset
        blocks = subset.blocks
}
