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
 *       one block, or against b_0) runs the full-pool perm-replay with export of the
 *       perm-discovery tables, then the ASR replay of the same genes
 *       (disambiguation_perms_main.py --detail-only), which gives the per-gene
 *       perm_pos_detail shards of the block; the b_0 tile also scores the b_0 labeling as
 *       full records (observed_b0_main.py) into b0_observed/.
 *    3. CAAS_CORE_OBSERVED writes the observed contract files (discovery.tab,
 *       background.output, background_genes.output, meta_caas/,
 *       ct_disambiguation/caas_convergence_master.csv) from the b_0 tiles' slices.
 *    4. CAAS_CORE_MERGE unions the tiles' shards (a gene's shards of the different blocks are
 *       joined into one, in the order of an unblocked replay) and derives the genome-wide
 *       tables and caas_perms.rds from them. It needs the gene universe, which the
 *       post-processing of the observed results provides, so it runs after
 *       CAAS_CORE_OBSERVED.
 *
 *  The aggregate RDS feeds the FCS p.perm path (fcs_enrich.R), which gives the CAAS
 *  scoring FCS report a permulation-corrected p.perm, as RER does.
 *
 *  Consumes:  resample directory (RESAMPLE), observed design (trait file or
 *             traitfile_H*.tab directory), alignments, species tree
 *  Produces:  caas_permulation/ (resample_perms.tab, caas_perms.rds,
 *             gene_cycle_scores.tsv, b0/, ...) and the observed contract files
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Helpers ──────────────────────────────────────────────────────────────────

// True unless params.caas_perms_postproc is off (Boolean false, or the strings false / 0 / no).
def caasPostprocOn() {
    def raw = params.containsKey('caas_perms_postproc') ? params.caas_perms_postproc : true
    return (raw instanceof Boolean) ? raw : !(raw?.toString()?.toLowerCase() in ['false', '0', 'no'])
}

// True when params.remove_intermediates is on: tasks delete the bulky intermediates that no later process reads.
def removeIntermediatesOn() {
    def raw = params.containsKey('remove_intermediates') ? params.remove_intermediates : false
    return (raw instanceof Boolean) ? raw : (raw?.toString()?.toLowerCase() in ['true', '1', 'yes'])
}

// True when the core also writes the strict tables: params.strict_conservation_output, and a discovery that tolerates
// conserved pairs (min_divergent_fraction < 1); with fraction 1 the strict rows are all the rows.
def strictOutputOn() {
    return (params.strict_conservation_output as boolean) && ((params.min_divergent_fraction ?: 0.5) as double) < 1.0
}

// CLI arguments that give disambiguation_perms_main.py the observed CT_POSTPROC filters (cluster trains and
// gene removal) for the per-cycle null pool. Empty when the filters are off or the gene annotation file is a
// NO_* sentinel. With params.caas_map_dir the trains measure their span in untrimmed columns, as the observed
// CT_FILTER does.
def caasPostprocArgs(gene_lengths) {
    if (!caasPostprocOn() || gene_lengths.name =~ /^NO_/) return ""
    def map_arg = params.caas_map_dir ? "--train-map-dir ${params.caas_map_dir}" : ""
    return "--postproc-filter --gene-lengths ${gene_lengths} --clust-minlen ${params.filter_minlen} --clust-maxcaas ${params.filter_maxcaas} --gene-filter-mode ${params.gene_filter_mode} --iqr-multiplier ${params.iqr_multiplier} --extreme-percentile ${params.extreme_threshold} ${params.remove_caas_clusters ? '' : '--keep-clusters'} ${map_arg}"
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
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'resample_perms.tab'

    input:
    path resample_dir
    val  n_perms
    val  block_size
    path caas_config   // observed design: trait file or multi-hypothesis dir (source of b_0)

    output:
    path "resample_perms.tab", emit: subset
    path "blocks/*", emit: blocks   // blocks/b0 and one blocks/blk_NNNN per block of cycles

    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    // A multi-hypothesis design (a directory of hypotheses) keeps all of them for b_0 even when the harvest holds no permuted cycle
    // (N = 0) and so wrote no fop_labelings.tab.
    def fop_design = "${params.multi_hypothesis}" == 'true'
    """
    if [ "${block_size}" -lt 1 ]; then echo "ERROR: caas_cycle_block_size must be at least 1 (got ${block_size})" >&2; exit 1; fi
    mkdir -p blocks/b0

    # Writes the cycles of \$1 (labelings in cycle order, tags b_N or b_N~H<m>) into blocks/blk_NNNN/resample_perms.tab and,
    # when \$2 is a pair table (header + rows keyed by base cycle), their rows into blocks/blk_NNNN/fop_pairs.tsv.
    split_blocks() {
        cut -f1 "\$1" | sed 's/~.*//' | awk '!seen[\$0]++' > block_cycles.txt
        n_cycles=\$(wc -l < block_cycles.txt)
        for (( k = 1; k <= (n_cycles + ${block_size} - 1) / ${block_size}; k++ )); do mkdir -p "\$(printf 'blocks/blk_%04d' \$k)"; done
        awk -F'\\t' -v B=${block_size} 'NR==FNR{blk[\$1]=int((FNR-1)/B)+1; next} {b=\$1; sub(/~.*/,"",b); print blk[b] "\\t" \$0}' block_cycles.txt "\$1" \\
            | awk -F'\\t' '{ if (\$1 != cur) { if (cur != "") close(f); cur = \$1; f = sprintf("blocks/blk_%04d/resample_perms.tab", cur) } sub(/^[0-9]+\\t/, ""); print > f }'
        if [ -n "\$2" ]; then
            hdr=\$(head -n 1 "\$2")
            awk -F'\\t' -v B=${block_size} 'NR==FNR{blk[\$1]=int((FNR-1)/B)+1; next} FNR>1 && (\$1 in blk){print blk[\$1] "\\t" \$0}' block_cycles.txt "\$2" \\
                | LC_ALL=C sort -s -k1,1n \\
                | awk -F'\\t' -v hdr="\$hdr" '{ if (\$1 != cur) { if (cur != "") close(f); cur = \$1; f = sprintf("blocks/blk_%04d/fop_pairs.tsv", cur); print hdr > f } sub(/^[0-9]+\\t/, ""); print > f }'
        fi
        rm -f block_cycles.txt
    }

    # Deterministic collect: cycles are ordered by their numeric id and the first N are kept
    # (N = 0 keeps none), so a resample with more cycles than N always yields the same null.
    # (-L: the resample directory is staged as a symlink.) The real labeling b_0 is rebuilt
    # from the observed design and prepended, so it is replayed through the same path as the
    # permulated cycles and split off downstream.
    FOP_TAB=""
    if [ -d "${resample_dir}" ]; then
        FOP_TAB=\$(find -L ${resample_dir} -name 'fop_labelings.tab' | head -n 1)
    fi
    # numeric cycle id of a row's tag: "b_12~H3" -> 12
    CYC='{b=\$1; sub(/^b_/,"",b); sub(/~.*/,"",b); print b"\\t"\$0}'

    if [ -n "\$FOP_TAB" ] || { [ "${fop_design}" = true ] && [ -d "${caas_config}" ]; }; then
        # FOP mirror: labelings are "<base>~H<m>". All hypothesis rows of the first N base
        # cycles are kept, with the matching fop_pairs.tsv rows (PSS weights for pooling).
        if [ -n "\$FOP_TAB" ]; then
            FOP_PAIRS=\$(find -L ${resample_dir} -name 'fop_pairs.tsv' | head -n 1)
            if [ -z "\$FOP_PAIRS" ]; then echo "ERROR: FOP labelings without fop_pairs.tsv (b_0 needs its PSS weights)" >&2; exit 1; fi
            awk -F'\\t' 'NF>=3 && \$1!~/^b_0(~|\$)/' "\$FOP_TAB" > candidates.tab
        else
            # no permuted cycle was harvested (N = 0): b_0 comes from the design alone
            : > candidates.tab
            printf 'cycle\\thypothesis_id\\tpair\\tspecies1\\tspecies2\\tpss_score\\n' > fop_pairs_header.tsv
            FOP_PAIRS=fop_pairs_header.tsv
        fi
        awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' candidates.tab | sort -u \\
            | awk '{b=\$1; sub(/^b_/,"",b); print b"\\t"\$1}' | sort -k1,1n | cut -f2 > base_all.txt
        awk -v n=${n_perms} 'n>0{print; if(++c>=n) exit}' base_all.txt > keep_base.txt
        awk -F'\\t' 'NR==FNR{k[\$1]=1; next} {b=\$1; sub(/~.*/,"",b); if(b in k) print}' \\
            keep_base.txt candidates.tab | awk -F'\\t' "\$CYC" | sort -s -k1,1n | cut -f2- > resample_perms.tab
        n=\$(awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' resample_perms.tab | sort -u | wc -l)
        rows=\$(wc -l < resample_perms.tab)
        avail=\$(wc -l < base_all.txt)
        echo "[caas_perms] FOP mirror: first \$n of \$avail base cycles in cycle order (\$rows hypothesis labelings; requested ${n_perms})"
        if [ "${n_perms}" -gt 0 ] && [ "\$rows" -eq 0 ]; then echo "ERROR: no FOP labelings selected" >&2; exit 1; fi
        split_blocks resample_perms.tab "\$FOP_PAIRS"
        ${run} $baseDir/subworkflows/CT/local/scripts/build_b0_labelings.py --config ${caas_config} --fop \\
            --labelings-out b0_labelings.tab --pairs-out b0_pairs.tsv
        cp b0_labelings.tab blocks/b0/resample_perms.tab
        { head -n 1 "\$FOP_PAIRS"; cat b0_pairs.tsv; } > blocks/b0/fop_pairs.tsv
        cat b0_labelings.tab resample_perms.tab > resample_perms.tmp && mv resample_perms.tmp resample_perms.tab
        # The selection files are not outputs; they run to hundreds of MB in FOP mode.
        rm -f all_resamples.tab candidates.tab base_all.txt keep_base.txt b0_labelings.tab b0_pairs.tsv fop_pairs_header.tsv
        exit 0
    fi

    if [ -d "${resample_dir}" ]; then
        cat \$(find -L ${resample_dir} -name 'resample_*.tab' | sort) > all_resamples.tab
    elif [ -f "${resample_dir}" ]; then
        cat ${resample_dir} > all_resamples.tab
    else
        echo "ERROR: resample input not found: ${resample_dir}" >&2; exit 1
    fi
    awk -F'\\t' 'NF>=3 && \$1!="b_0"' all_resamples.tab | awk -F'\\t' "\$CYC" | sort -s -k1,1n | cut -f2- > candidates.tab
    awk -v n=${n_perms} 'n>0{print; if(++c>=n) exit}' candidates.tab > resample_perms.tab
    n=\$(wc -l < resample_perms.tab)
    avail=\$(wc -l < candidates.tab)
    echo "[caas_perms] first \$n of \$avail permuted labelings in cycle order (requested ${n_perms})"
    if [ "${n_perms}" -gt 0 ] && [ "\$n" -eq 0 ]; then echo "ERROR: no permuted labelings selected" >&2; exit 1; fi
    split_blocks resample_perms.tab ""
    ${run} $baseDir/subworkflows/CT/local/scripts/build_b0_labelings.py --config ${caas_config} --labelings-out b0_labelings.tab
    cp b0_labelings.tab blocks/b0/resample_perms.tab
    cat b0_labelings.tab resample_perms.tab > resample_perms.tmp && mv resample_perms.tmp resample_perms.tab
    rm -f all_resamples.tab candidates.tab b0_labelings.tab
    """
}

// ── 2. One task per (gene batch, cycle block): perm-replay, then the ASR replay of the same genes ───
// A tile is the replay of one batch of genes against the labelings of one cycle block (or of b_0, the tile
// whose block directory is blocks/b0). The task replays the batch's alignments through `ct perm-replay` (full
// position pool, perm-discovery export; each ct process parses the block's labelings once for all the genes it
// takes) and feeds the exports straight into disambiguation_perms_main.py --detail-only, so no gene waits for the
// rest of the genome between the two steps. Pass B scores each cycle against genome-wide pools and runs once, in
// CAAS_CORE_MERGE, over the union of the tiles' shards.
// The b_0 tile also exports the b_0 slice of each gene (the discovery.tab rows of the real labeling and the
// positions it tested); observed_b0_main.py then scores those rows as full records.
// Output: perm_pos_detail/ (one gz shard per gene with a detection in the block; the b_0 shards under b0/) and, for the
// b_0 tile, b0_observed/ (per gene with b_0 hits: its discovery rows, its tested positions and its master rows).
process CAAS_CORE_BATCHED {
    tag "$batchID (${batchSize} genes) ${block.name}"
    label 'process_resample'

    input:
    tuple val(batchID), val(batchSize), val(manifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(block)   // block: blocks/blk_NNNN or blocks/b0 of SUBSET_RESAMPLE_PERMS
    file caas_config
    path tree_file
    path gene_lengths // gene_ensembl_file (CT_POSTPROC filters) or NO_FILE
    path taxid_map    // tax_id map of the species (curated by NAME_CURATION, or params.tax_id) or the NO_FILE sentinel

    output:
    path "perm_pos_detail", emit: pos_detail
    path "b0_observed", emit: b0_observed, optional: true                    // b_0 tile only
    path "perm_pos_detail_strict", emit: pos_detail_strict, optional: true   // strict pooling of the same replay (strictOutputOn)
    path "b0_observed_strict", emit: b0_observed_strict, optional: true      // b_0 rows without conserved pairs, scored on their own (b_0 tile only)

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct"
        : "$baseDir/subworkflows/CT/local/ct"
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = taxid_map.name != 'NO_FILE' ? taxid_map : ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def max_tasks_per_child = params.ct_disambig_max_tasks_per_child ?: 50
    def postproc_args = caasPostprocArgs(gene_lengths)
    // The perm-discovery exports (~13-16 MB per gene at K=3) are read only by the ASR replay below; with
    // remove_intermediates they are deleted right after it.
    def prune_exports = removeIntermediatesOn()
    // With remove_intermediates the shards are moved out of caas_perms_out instead of duplicated.
    def shard_move = removeIntermediatesOn()
    // The strict tables come from the same replay: the null shards pooled without the hypotheses that have conserved
    // pairs, and the b_0 rows without conserved pairs scored again (the observed labeling is a single labeling).
    def strict_on = strictOutputOn()
    // The b_0 tile replays the real labeling alone; the block tiles hold permuted cycles only.
    def is_b0 = (block.name == 'b0')
    def strict_flag = strict_on ? '--strict-output' : ''
    def strict_observed_flag = strict_on ? '--strict-output-dir b0_observed_strict' : ''
    def strict_empty = strict_on ? 'mkdir -p caas_perms_out/perm_pos_detail_strict' : ''
    // The observed labeling as full records: one master shard per gene with b_0 hits (the b_0 tile only).
    def observed_b0 = !is_b0 ? '' : """
    # The observed labeling as full records: one master shard per gene with b_0 hits.
    if compgen -G 'b0_observed/*.b0.discovery.tsv' > /dev/null; then
        ${run} ./observed_b0_main.py \\
            --alignment-dir ${params.alignment} \\
            --tree ${tree_file} \\
            --b0-dir b0_observed ${strict_observed_flag} \\
            --design ${caas_config} \\
            \$FOP_PAIRS_ARG \\
            --output-dir b0_observed \\
            --asr-model ${params.ct_disambig_asr_model} \\
            --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
            --workers ${task.cpus} \\
            --max-tasks-per-child ${max_tasks_per_child} \\
            --asr-cache-dir ${asr_cache_dir} \\
            ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
            ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    fi"""
    def strict_move = strict_on ? """
    ${shard_move ? 'mv' : 'cp -R'} caas_perms_out/perm_pos_detail_strict perm_pos_detail_strict
    if [ -d caas_perms_out/b0/perm_pos_detail_strict ]; then ${shard_move ? 'mv' : 'cp -R'} caas_perms_out/b0/perm_pos_detail_strict perm_pos_detail_strict/b0; fi""" : ''
    """
    # One worker pool per step runs the genes in parallel; BLAS and OpenMP stay at one thread each.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    # Memory of one replay worker: the ASR replay splits a gene whose replay would not fit in it into cycle chunks.
    export CAAS_PERMS_WORKER_MEM_MB=\$(( ${task.memory.toMega()} / ${task.cpus} ))
    mkdir -p perm_disc ${is_b0 ? 'b0_observed' : ''}
    ${is_b0 && strict_on ? 'mkdir -p b0_observed_strict' : ''}
    # Bounds the chunk buffers of the vectorized kernel, per call, to a quarter of one worker's share of task.memory.
    export CT_PERM_REPLAY_CHUNK_MEM_MB=\$(( ${task.memory.toMega()} / ${task.cpus} / 4 ))
    cat > ${batchID}.manifest.tsv <<'EOF'
""" + manifestText + """EOF

    # With multi_hypothesis the design is a directory of traitfile_H*.tab (same K): one .tab gives the pair count
    if [ -d "${caas_config}" ]; then
        _cfg_file=\$(find -L ${caas_config} -type f -name '*.tab' | sort | head -n 1)
    else
        _cfg_file="${caas_config}"
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

    bash $baseDir/subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh \\
        --batch-id ${batchID} \\
        --manifest ${batchID}.manifest.tsv \\
        --caas-config ${caas_config} \\
        --resampled-path ${block}/resample_perms.tab \\
        --workers ${task.cpus} \\
        --ali-format ${params.ali_format} \\
        --ct-bin ${ctBinary} \\
        --export-groups 0 \\
        --export-perm-discovery 1 \\
        --export-b0 ${is_b0 ? 1 : 0} \\
        --extra-args-file .ct_perm_replay_batch_args
    find . -maxdepth 1 -name '*.perm_replay.discovery.output' -exec mv {} perm_disc/ \\;
    ${is_b0 ? "find . -maxdepth 1 \\( -name '*.b0.discovery.tsv' -o -name '*.b0.background' \\) -exec mv {} b0_observed/ \\;" : ''}

    # The PSS weights of the block's cycles (absent for a design without FOP hypotheses).
    FOP_PAIRS_ARG=""
    if [ -f ${block}/fop_pairs.tsv ]; then FOP_PAIRS_ARG="--fop-pairs ${block}/fop_pairs.tsv"; fi

    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    mkdir -p caas_perms_out
    # A batch whose genes have no CAAS in any cycle exports nothing: its shard directory is empty.
    if [ -n "\$(ls -A perm_disc)" ]; then
        ${run} ./disambiguation_perms_main.py \\
            --alignment-dir ${params.alignment} \\
            --tree ${tree_file} \\
            --perm-discovery perm_disc \\
            --resample-dir ${block} \\
            --output-dir caas_perms_out \\
            --detail-only ${strict_flag} \\
            \$FOP_PAIRS_ARG ${postproc_args} \\
            --asr-model ${params.ct_disambig_asr_model} \\
            --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
            --workers ${task.cpus} \\
            --max-tasks-per-child ${max_tasks_per_child} \\
            --asr-cache-dir ${asr_cache_dir} \\
            --seed ${params.seed ?: 1998} \\
            ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
            ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    else
        mkdir -p caas_perms_out/perm_pos_detail
        ${strict_empty}
    fi
    ${prune_exports ? 'rm -rf perm_disc/*.perm_replay.discovery.output' : ''}
    ${shard_move ? 'mv' : 'cp -R'} caas_perms_out/perm_pos_detail perm_pos_detail
    # b_0 shards sit in a subdirectory the null readers never glob; CAAS_CORE_MERGE scores them as a one-labeling run.
    if [ -d caas_perms_out/b0/perm_pos_detail ]; then ${shard_move ? 'mv' : 'cp -R'} caas_perms_out/b0/perm_pos_detail perm_pos_detail/b0; fi
    ${strict_move}
    ${shard_move ? 'rm -rf caas_perms_out' : ''}

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
        // Every gene batch replays against every block directory.
        def tiles = live.combine(blocks.flatten())
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
        pos_detail     = core.pos_detail
        b0_observed    = core.b0_observed
        pos_detail_strict = core.pos_detail_strict   // empty unless the strict tables are on
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
    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    ${py} ./contract_main.py --b0-dirs b0obs_* --design ${design} --output-dir .
    """
}

// ── 4. Merge the tiles' shards and derive the genome-wide null (no ASR replay) ───
// The tiles' shard directories are united with hard links (copies where the filesystem
// refuses them) into perm_pos_detail/, which CT_ACCUMULATION also reads; a gene replayed in
// several cycle blocks has its shards joined into one, in the cycle order of an unblocked
// replay. A single perm_pos_detail.tsv.gz instead of directories is read as it is.
//
// Each gene's score is calibrated against the genome-wide pool of position scores of its
// cycle, so pass B runs once over the union, through reaggregate_perm_scores.py (the
// aggregation of gene_wrapper.py itself, not a second implementation: the gene statistic is
// F(max)^n over heavily tied values, so a 1e-16 difference in how the per-position sum
// accumulates can flip a tie and the power n amplifies it). Pass B costs minutes, against
// the hours of the ASR replay, so caas_perms.rds is always rebuilt from the shards: the null
// must hold the same gene statistic as the observed gene score, or the FCS p.perm would
// compare two different quantities.
process CAAS_CORE_MERGE {
    tag "caas_core_merge"
    label 'process_medium'
    // ext.publish_prefix is "/strict" for the strict alias (CAAS_CORE_MERGE_STRICT)
    publishDir path: { "${params.outdir}${task.ext.publish_prefix ?: ''}/caas_permulation" }, mode: 'copy', overwrite: true,
               pattern: '{caas_perms.rds,gene_cycle_scores.tsv,perm_pos_sample.tsv,perm_pos_quantiles.tsv,perm_pos_cycle_caas.tsv.gz,cycle_roster.txt,removed_units.tsv,b0}'

    input:
    path batchDetail, stageAs: 'batch_*'   // the tiles' shard directories, or one legacy perm_pos_detail.tsv.gz
    path universe
    path gene_lengths   // gene_ensembl_file (gene removal) or NO_FILE
    path labelings      // the labelings file the cycles were replayed from, or NO_FILE: N is then read from the detail rows

    output:
    path "perm_pos_detail",           emit: pos_detail, optional: true   // the union of the batches' shards (absent for a legacy file)
    path "b0",                        emit: b0, optional: true   // the real labeling rebuilt like the null
    path "caas_perms.rds",            emit: perms
    path "gene_cycle_scores.tsv",     emit: gene_cycle_scores
    path "removed_units.tsv",         emit: removed_units, optional: true   // null gene removal (only when a gene annotation was given)
    path "cycle_roster.txt",          emit: cycle_roster,  optional: true   // every replayed cycle, the ones with no detection included: N of every consumer of the null (absent when the labelings file is NO_LABELINGS)
    path "perm_pos_cycle_caas.tsv.gz", emit: pos_cycle_caas, optional: true
    path "perm_pos_sample.tsv",       emit: pos_sample,    optional: true
    path "perm_pos_quantiles.tsv",    emit: pos_quantiles, optional: true

    script:
    def disambig_local = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scoring_local  = "${baseDir}/subworkflows/SCORING/local"
    // The sentinel convention is a NO_* basename (NO_FILE here, NO_BACKGROUND from SCORING's
    // resolved_background), so the prefix is matched, not one literal.
    def universe_arg   = universe.name.startsWith('NO_') ? "" : "--universe ${universe}"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript'  : 'Rscript'
    // N = 0 (b_0 only): the null has no shard, and says so explicitly. Written without `?:`, because Groovy reads a numeric 0 as false.
    def empty_null_arg = (params.caas_full_perms != null && (params.caas_full_perms as int) == 0) ? "--empty-null" : ""
    // N is the number of cycles replayed, which the detail rows alone understate when a cycle left no row.
    // Gene removal needs the annotation file, so it is on only when gene_lengths is not a NO_* sentinel.
    def roster_arg = labelings.name.startsWith('NO_') ? "" : "--cycles-from ${labelings}"
    def map_arg = params.caas_map_dir ? "--map-dir '${params.caas_map_dir}'" : ""
    def removal_args = (caasPostprocOn() && !gene_lengths.name.startsWith('NO_')) ? "--gene-lengths ${gene_lengths} ${map_arg} --gene-filter-mode ${params.gene_filter_mode} --iqr-multiplier ${params.iqr_multiplier} --extreme-percentile ${params.extreme_threshold} ${params.remove_caas_clusters ? '' : '--keep-clusters'}" : ""
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    # reaggregate_perm_scores.py imports src.utils.gene_wrapper relative to its own location.
    cp -R ${disambig_local}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    entries=(batch_*)
    if [ "\${#entries[@]}" -eq 1 ] && [ -f "\${entries[0]}" ]; then
        DETAIL="\${entries[0]}"
    else
        # Hard links when the first shard can be linked from here, copies otherwise; the shards of a gene
        # replayed in several cycle blocks are joined in the cycle order of an unblocked replay (union_shards.sh).
        UNION_PYTHON="${py}" bash $baseDir/subworkflows/CT/local/scripts/union_shards.sh perm_pos_detail "\${entries[@]}"
        DETAIL=perm_pos_detail
    fi

    ${py} ./reaggregate_perm_scores.py \\
        --detail "\$DETAIL" \\
        --output-dir . \\
        --seed ${params.seed ?: 1998} ${empty_null_arg} ${roster_arg} ${removal_args}

    # b_0 rebuilt from its merged shards as a one-labeling run: same code, own rank/size pools.
    if [ -d "\$DETAIL/b0" ] && compgen -G "\$DETAIL/b0/*.tsv.gz" > /dev/null; then   # no shard: nothing was detected for the real labeling
        mkdir -p b0
        ${py} ./reaggregate_perm_scores.py \\
            --detail "\$DETAIL/b0" \\
            --output-dir b0 \\
            --seed ${params.seed ?: 1998} ${removal_args}
        # The b_0 per-gene shards stay next to its scores
        cp -Rl "\$DETAIL/b0" b0/perm_pos_detail 2>/dev/null || cp -RL "\$DETAIL/b0" b0/perm_pos_detail
    fi

    ROSTER_ARG=""
    if [ -f cycle_roster.txt ]; then ROSTER_ARG="--cycles cycle_roster.txt"; fi
    ${rs} ${scoring_local}/src/scoring_caas_perms.R \\
        --gene-cycle-scores gene_cycle_scores.tsv \\
        \$ROSTER_ARG \\
        ${universe_arg} \\
        --output caas_perms.rds

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
