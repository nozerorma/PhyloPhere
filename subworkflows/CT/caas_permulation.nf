#!/usr/bin/env nextflow

/*
 * CAAS permulation-excess null
 * ────────────────────────────
 * Builds a genome-wide *excess* null for CAAS FCS pathway enrichment:
 *   1. SUBSET_RESAMPLE_PERMS: take the first N (caas_full_perms) permuted
 *      labelings in cycle order (N = 0 keeps none). The real labeling b_0 is rebuilt from the
 *      observed design and replayed alongside; it is never part of the null (its shards go to
 *      caas_permulation/b0/).
 *   2. CAAS_CORE_BATCHED (one task per batch of ct_core_batch_size genes): full-pool
 *      perm-replay with export_perm_discovery, then the ASR replay of the same genes
 *      (disambiguation_perms_main.py --detail-only) → per-gene perm_pos_detail shards, and the
 *      observed labeling (b_0) as full records (observed_b0_main.py) → b0_observed/.
 *      Perm-discovery exports that already exist are disambiguated without a replay.
 *   3. CAAS_CORE_OBSERVED writes the observed contract files (discovery.tab, background.output,
 *      background_genes.output, meta_caas/, ct_disambiguation/caas_convergence_master.csv) from the
 *      batches' b_0 slices.
 *   4. CAAS_CORE_MERGE unions the batches' shards and derives the genome-wide tables and
 *      caas_perms.rds from them. It needs the gene universe, which the post-processing of the
 *      observed results provides, so it runs after CAAS_CORE_OBSERVED.
 *
 * The aggregate RDS feeds the existing FCS p.perm path (fcs_enrich.R), giving
 * the CAAS scoring FCS report a permulation-corrected p.perm — exactly like RER.
 * See docs/CAAS_PERMULATION_EXCESS.md.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

// True unless params.caas_perms_postproc is off (Boolean false, or the strings false / 0 / no).
def caasPostprocOn() {
    def raw = params.containsKey('caas_perms_postproc') ? params.caas_perms_postproc : true
    return (raw instanceof Boolean) ? raw : !(raw?.toString()?.toLowerCase() in ['false', '0', 'no'])
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

// ── 1. Collect the first N permuted labelings in cycle order (+ optional b_0) ─
process SUBSET_RESAMPLE_PERMS {
    tag "caas_perms_subset|N=${n_perms}"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'resample_perms.tab'

    input:
    path resample_dir
    val  n_perms
    path caas_config   // observed design: trait file or multi-hypothesis dir (source of b_0)

    output:
    path "resample_perms.tab", emit: subset
    path "fop_pairs.tsv", emit: fop_pairs, optional: true

    script:
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    # Deterministic collect: cycles are ordered by their numeric id and the first N are
    # kept (N = 0 keeps none), so a resample reused with more cycles than N yields the same null
    # every time. (-L: the resample dir is staged as a symlink.) The real labeling b_0 is rebuilt
    # from the observed design and prepended, so it is replayed through the same path as the
    # permuted cycles and split off downstream.
    FOP_TAB=""
    if [ -d "${resample_dir}" ]; then
        FOP_TAB=\$(find -L ${resample_dir} -name 'fop_labelings.tab' | head -n 1)
    fi
    # numeric cycle id of a row's tag: "b_12~H3" -> 12
    CYC='{b=\$1; sub(/^b_/,"",b); sub(/~.*/,"",b); print b"\\t"\$0}'

    if [ -n "\$FOP_TAB" ]; then
        # FOP mirror: labelings are "<base>~H<m>". Keep ALL hypothesis rows of the first N
        # base cycles, plus the matching fop_pairs.tsv rows (PSS weights for pooling).
        FOP_PAIRS=\$(find -L ${resample_dir} -name 'fop_pairs.tsv' | head -n 1)
        if [ -z "\$FOP_PAIRS" ]; then echo "ERROR: FOP labelings without fop_pairs.tsv (b_0 needs its PSS weights)" >&2; exit 1; fi
        awk -F'\\t' 'NF>=3 && \$1!~/^b_0(~|\$)/' "\$FOP_TAB" > candidates.tab
        awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' candidates.tab | sort -u \\
            | awk '{b=\$1; sub(/^b_/,"",b); print b"\\t"\$1}' | sort -k1,1n | cut -f2 > base_all.txt
        awk -v n=${n_perms} 'n>0{print; if(++c>=n) exit}' base_all.txt > keep_base.txt
        awk -F'\\t' 'NR==FNR{k[\$1]=1; next} {b=\$1; sub(/~.*/,"",b); if(b in k) print}' \\
            keep_base.txt candidates.tab | awk -F'\\t' "\$CYC" | sort -s -k1,1n | cut -f2- > resample_perms.tab
        head -n 1 "\$FOP_PAIRS" > fop_pairs.tsv
        awk -F'\\t' 'NR==FNR{k[\$1]=1; next} FNR>1 && (\$1 in k){print}' keep_base.txt "\$FOP_PAIRS" >> fop_pairs.tsv
        n=\$(awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' resample_perms.tab | sort -u | wc -l)
        rows=\$(wc -l < resample_perms.tab)
        avail=\$(wc -l < base_all.txt)
        echo "[caas_perms] FOP mirror: first \$n of \$avail base cycles in cycle order (\$rows hypothesis labelings; requested ${n_perms})"
        if [ "${n_perms}" -gt 0 ] && [ "\$rows" -eq 0 ]; then echo "ERROR: no FOP labelings selected" >&2; exit 1; fi
        ${run} $baseDir/subworkflows/CT/local/scripts/build_b0_labelings.py --config ${caas_config} --fop \\
            --labelings-out b0_labelings.tab --pairs-out b0_pairs.tsv
        cat b0_labelings.tab resample_perms.tab > resample_perms.tmp && mv resample_perms.tmp resample_perms.tab
        cat b0_pairs.tsv >> fop_pairs.tsv
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
    ${run} $baseDir/subworkflows/CT/local/scripts/build_b0_labelings.py --config ${caas_config} --labelings-out b0_labelings.tab
    cat b0_labelings.tab resample_perms.tab > resample_perms.tmp && mv resample_perms.tmp resample_perms.tab
    """
}

// ── 2. One task per gene batch: perm-replay, then the ASR replay of the same genes ─────────────────────
// The task replays the batch's alignments through `ct perm-replay` (full position pool, perm-discovery
// export) and feeds the exports straight into disambiguation_perms_main.py --detail-only, so no gene waits
// for the rest of the genome between the two steps. Pass B scores each cycle against genome-wide pools
// and runs once, in CAAS_CORE_MERGE, over the union of the batches' shards.
// A batch that arrives with perm-discovery files already computed (permDiscFiles) skips the replay.
// The same replay exports the b_0 slice of each gene (the discovery.tab rows of the real labeling and the
// positions it tested); observed_b0_main.py then scores those rows as full records.
// Output: perm_pos_detail/ (one gz shard per gene; the b_0 shards under b0/), b0_observed/ (per gene with b_0
// hits: its discovery rows, its tested positions and its master rows), and the perm-discovery exports of a
// replayed batch. A batch that reuses exports has no b_0 slice: b0_observed/ is empty.
process CAAS_CORE_BATCHED {
    tag "$batchID (${batchSize} genes)"
    label 'process_resample'
    publishDir path: "${params.outdir}/caas_permulation/perm_disc", mode: 'copy', overwrite: true,
               pattern: 'perm_disc/*.perm_replay.discovery.output', saveAs: { fn -> fn.tokenize('/').last() }

    input:
    tuple val(batchID), val(batchSize), val(manifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(permDiscFiles, stageAs: 'perm_disc_in/*')
    path resample_subset
    file caas_config
    path tree_file
    path fop_pairs    // fop_pairs.tsv (FOP mirror) or NO_FOP_PAIRS sentinel
    path gene_lengths // gene_ensembl_file (CT_POSTPROC filters) or NO_FILE

    output:
    path "perm_pos_detail", emit: pos_detail
    path "b0_observed", emit: b0_observed
    path "perm_disc/*.perm_replay.discovery.output", emit: perm_discovery, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct"
        : "$baseDir/subworkflows/CT/local/ct"
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def replay = manifestText.trim() ? true : false
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def max_tasks_per_child = params.ct_disambig_max_tasks_per_child ?: 50
    def postproc_args = caasPostprocArgs(gene_lengths)
    """
    # One worker pool per step runs the genes in parallel; BLAS and OpenMP stay at one thread each.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    mkdir -p perm_disc b0_observed
    PERM_DISC=perm_disc_in
    if ${replay}; then
        # Bounds the vectorized kernel's per-call chunk buffers to a quarter of one worker's share of task.memory.
        export CT_PERM_REPLAY_CHUNK_MEM_MB=\$(( ${task.memory.toMega()} / ${task.cpus} / 4 ))
        cat > ${batchID}.manifest.tsv <<'EOF'
""" + manifestText + """EOF

        # multi_hypothesis mode passes a directory of traitfile_H*.tab (same K); resolve it to one .tab
        if [ -d "${caas_config}" ]; then
            _cfg_file=\$(find -L ${caas_config} -type f -name '*.tab' | sort | head -n 1)
        else
            _cfg_file="${caas_config}"
        fi
        n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$_cfg_file" | sort -nu | wc -l | tr -d ' ')
        _frac() { awk -v n="\$n_pairs" -v f="\$1" 'BEGIN{printf "%d", int(n*f)}'; }
        declare -a extra_opts=(--patterns "${params.patterns}")
        if [ "${params.miss_pair}" = "true" ]; then extra_opts+=(--miss_pair); fi
        if [ "${params.caap_mode}" = "true" ]; then extra_opts+=(--caap_mode); fi
        extra_opts+=(--max_conserved \$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}'))
        extra_opts+=(--max_bg_gaps \$(_frac ${params.max_bg_gaps_fraction}) --max_fg_gaps \$(_frac ${params.max_fg_gaps_fraction}) --max_gaps \$(_frac ${params.max_gaps_fraction}))
        extra_opts+=(--max_bg_miss \$(_frac ${params.max_bg_miss_fraction}) --max_fg_miss \$(_frac ${params.max_fg_miss_fraction}) --max_miss \$(_frac ${params.max_miss_fraction}))
        echo "\${extra_opts[@]}" > .ct_perm_replay_batch_args

        bash $baseDir/subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh \\
            --batch-id ${batchID} \\
            --manifest ${batchID}.manifest.tsv \\
            --caas-config ${caas_config} \\
            --resampled-path ${resample_subset} \\
            --workers ${task.cpus} \\
            --ali-format ${params.ali_format} \\
            --ct-bin ${ctBinary} \\
            --progress-log 0 \\
            --export-groups 0 \\
            --export-perm-discovery 1 \\
            --export-b0 1 \\
            --extra-args-file .ct_perm_replay_batch_args
        find . -maxdepth 1 -name '*.perm_replay.discovery.output' -exec mv {} perm_disc/ \\;
        find . -maxdepth 1 \\( -name '*.b0.discovery.tsv' -o -name '*.b0.background' \\) -exec mv {} b0_observed/ \\;
        PERM_DISC=perm_disc
    fi

    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    mkdir -p caas_perms_out
    # A batch whose genes have no CAAS in any cycle exports nothing: its shard directory is empty.
    if [ -n "\$(ls -A \$PERM_DISC)" ]; then
        ${run} ./disambiguation_perms_main.py \\
            --alignment-dir ${params.alignment} \\
            --tree ${tree_file} \\
            --perm-discovery \$PERM_DISC \\
            --resample-dir . \\
            --output-dir caas_perms_out \\
            --detail-only \\
            ${fop_pairs.name =~ /^NO_/ ? '' : "--fop-pairs ${fop_pairs}"} ${postproc_args} \\
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
    fi
    cp -R caas_perms_out/perm_pos_detail perm_pos_detail
    # b_0 shards sit in a subdirectory the null readers never glob; CAAS_CORE_MERGE scores them as a one-labeling run.
    if [ -d caas_perms_out/b0/perm_pos_detail ]; then cp -R caas_perms_out/b0/perm_pos_detail perm_pos_detail/b0; fi

    # The observed labeling as full records: one master shard per gene with b_0 hits.
    if compgen -G 'b0_observed/*.b0.discovery.tsv' > /dev/null; then
        ${run} ./observed_b0_main.py \\
            --alignment-dir ${params.alignment} \\
            --tree ${tree_file} \\
            --b0-dir b0_observed \\
            --design ${caas_config} \\
            ${fop_pairs.name =~ /^NO_/ ? '' : "--fop-pairs ${fop_pairs}"} \\
            --output-dir b0_observed \\
            --asr-model ${params.ct_disambig_asr_model} \\
            --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
            --workers ${task.cpus} \\
            --max-tasks-per-child ${max_tasks_per_child} \\
            --asr-cache-dir ${asr_cache_dir} \\
            ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
            ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    fi
    """
}

// Batches of `ct_core_batch_size` genes (1 = one task per gene), in gene-name order. Alignments to replay,
// or perm-discovery files already computed (reuse), arrive on separate channels; each batch carries its own
// kind and the other slot holds a sentinel.
workflow CAAS_CORE {
    take:
        align_tuple      // Channel<tuple(id, alignmentFile)> to replay (empty when reusing exports)
        reuse_disc       // Channel<List<file>> of perm-discovery exports to disambiguate (empty when replaying)
        caas_config
        resample_subset
        tree_file        // gated on the ASR cache when ASR is computed live
        fop_pairs
        gene_lengths

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
                tuple(id, batch.size(), manifest, batch.collect { row -> row[1] }.unique { f -> f.name }, file('NO_PERM_DISC'))
            }
        def reuseCounter = 0
        def reuse = reuse_disc
            .flatMap { files -> files.sort { f -> f.name } }
            .collate(batchSize)
            .map { batch ->
                def id = String.format('caas_core_reuse_batch_%05d', ++reuseCounter)
                tuple(id, batch.size(), '', file('NO_ALIGNMENTS'), batch)
            }
        // The single-item inputs (resample subset, FOP pairs, gene lengths, trait file, tree) are queue channels
        // once they have crossed a take:/emit: boundary, not Nextflow value channels. A process pairs its input
        // channels positionally and stops when any one runs out, so paired against the many-item batch channel
        // each of them would end CAAS_CORE_BATCHED after its first batch. .collect().map { items -> items[0] }
        // makes each a value channel that every batch can read.
        def subset_bc = resample_subset.collect().map { items -> items[0] }
        def fop_bc    = fop_pairs.collect().map { items -> items[0] }
        def lengths_bc = gene_lengths.collect().map { items -> items[0] }
        def config_bc = caas_config.collect().map { items -> items[0] }
        def tree_bc   = tree_file.collect().map { items -> items[0] }
        def core = CAAS_CORE_BATCHED(live.mix(reuse), subset_bc, config_bc, tree_bc, fop_bc, lengths_bc)

    emit:
        pos_detail     = core.pos_detail
        b0_observed    = core.b0_observed
        perm_discovery = core.perm_discovery
}

// ── 3. The observed contract files, from the b_0 slices of the batches ───────────────────────────────
// contract_main.py turns the b0_observed/ directories into discovery.tab, background.output,
// background_genes.output, meta_caas/ and ct_disambiguation/caas_convergence_master.csv: every table holds the
// genes in name order, whatever order the batches finished in. Batches that reuse perm-discovery exports carry
// no b_0 slice and give no file. Published where the rest of the pipeline reads them: caastools/, meta_caas/meta_caas/
// and ct_disambiguation/.
process CAAS_CORE_OBSERVED {
    tag "caas_core_observed"
    label 'process_low'
    publishDir path: "${params.outdir}/caastools", mode: 'copy', overwrite: true, pattern: '{discovery.tab,background.output,background_genes.output}'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: 'meta_caas/**'
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, pattern: 'ct_disambiguation/**'

    input:
    path b0Observed, stageAs: 'b0obs_*'   // the batches' b0_observed directories
    path design                           // observed design: the master columns come from it

    output:
    path "discovery.tab",                emit: discovery, optional: true
    path "background.output",            emit: background, optional: true
    path "background_genes.output",      emit: background_genes, optional: true
    path "ct_disambiguation",            emit: results_dir, optional: true
    path "ct_disambiguation/caas_convergence_master.csv", emit: master_csv, optional: true
    path "meta_caas",                    emit: meta_caas, optional: true
    path "meta_caas/global_meta_caas.tsv", emit: global_meta_caas, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    ${py} ./contract_main.py --b0-dirs b0obs_* --design ${design} --output-dir .
    """
}

// ── 3. Merge the batches' shards and derive the genome-wide null (no ASR replay) ──────────────────────
// Batches partition genes disjointly, so the shard directories are united with hard links (copies where the
// filesystem refuses them) into perm_pos_detail/, which CT_ACCUMULATION also reads. A single legacy
// perm_pos_detail.tsv.gz is read as it is.
//
// Each gene's score is calibrated against its cycle's genome-wide pool of position scores, so pass B runs once
// over the union, through reaggregate_perm_scores.py (gene_wrapper.py's own aggregation, not a second
// implementation: the gene statistic is F(max)^n over heavily tied values, so a 1e-16 difference in how the
// per-position sum accumulates can flip a tie and the ^n amplifies it). Pass B costs minutes against the hours
// of the ASR replay, which is why a caas_perms.rds imported from an earlier run is never trusted: it holds
// whatever gene statistic was current when it was built, and the null must hold the same statistic as the
// observed gene score or the FCS p.perm compares two different quantities.
process CAAS_CORE_MERGE {
    tag "caas_core_merge"
    label 'process_medium'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true,
               pattern: '{caas_perms.rds,gene_cycle_scores.tsv,perm_pos_sample.tsv,perm_pos_quantiles.tsv,perm_pos_cycle_caas.tsv.gz,removed_units.tsv,b0}'

    input:
    path batchDetail, stageAs: 'batch_*'   // batch shard directories, or one legacy perm_pos_detail.tsv.gz
    path universe
    path gene_lengths   // gene_ensembl_file (gene removal) or NO_FILE

    output:
    path "perm_pos_detail",           emit: pos_detail, optional: true   // the union of the batches' shards (absent for a legacy file)
    path "b0",                        emit: b0, optional: true   // the real labeling rebuilt like the null
    path "caas_perms.rds",            emit: perms
    path "gene_cycle_scores.tsv",     emit: gene_cycle_scores
    path "removed_units.tsv",         emit: removed_units, optional: true   // null gene removal (only when a gene annotation was given)
    path "perm_pos_cycle_caas.tsv.gz", emit: pos_cycle_caas, optional: true
    path "perm_pos_sample.tsv",       emit: pos_sample,    optional: true
    path "perm_pos_quantiles.tsv",    emit: pos_quantiles, optional: true

    script:
    def disambig_local = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scoring_local  = "${baseDir}/subworkflows/SCORING/local"
    // Sentinel convention is a NO_* basename (NO_FILE here, NO_BACKGROUND from
    // SCORING's resolved_background) — match the prefix, not one literal.
    def universe_arg   = universe.name.startsWith('NO_') ? "" : "--universe ${universe}"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript'  : 'Rscript'
    // Gene removal needs the annotation file.
    def removal_args = (caasPostprocOn() && !gene_lengths.name.startsWith('NO_')) ? "--gene-lengths ${gene_lengths} --gene-filter-mode ${params.gene_filter_mode} --iqr-multiplier ${params.iqr_multiplier} --extreme-percentile ${params.extreme_threshold} ${params.remove_caas_clusters ? '' : '--keep-clusters'}" : ""
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    # reaggregate_perm_scores.py imports src.utils.gene_wrapper relative to its own location.
    cp -R ${disambig_local}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    entries=(batch_*)
    if [ "\${#entries[@]}" -eq 1 ] && [ -f "\${entries[0]}" ]; then
        DETAIL="\${entries[0]}"
    else
        # Hard links when the first shard can be linked from here, copies otherwise.
        bash $baseDir/subworkflows/CT/local/scripts/union_shards.sh perm_pos_detail "\${entries[@]}"
        DETAIL=perm_pos_detail
    fi

    ${py} ./reaggregate_perm_scores.py \\
        --detail "\$DETAIL" \\
        --output-dir . \\
        --seed ${params.seed ?: 1998} ${removal_args}

    # b_0 rebuilt from its merged shards as a one-labeling run: same code, own rank/size pools.
    if [ -d "\$DETAIL/b0" ]; then
        mkdir -p b0
        ${py} ./reaggregate_perm_scores.py \\
            --detail "\$DETAIL/b0" \\
            --output-dir b0 \\
            --seed ${params.seed ?: 1998} ${removal_args}
        # keep b_0's own per-gene shards next to its scores: compare_b0.py reads them (checkpoints B, C)
        cp -RL "\$DETAIL/b0" b0/perm_pos_detail
    fi

    ${rs} ${scoring_local}/src/scoring_caas_perms.R \\
        --gene-cycle-scores gene_cycle_scores.tsv \\
        ${universe_arg} \\
        --output caas_perms.rds

    """
}

// ── Subworkflow: resample subset (runs INSIDE ct.nf, where the resample directory is available) ───────
workflow CAAS_PERMS_PREP {
    take:
        caas_config        // path (trait file)
        resample_dir       // path (resample_*.tab directory)

    main:
        def subset = SUBSET_RESAMPLE_PERMS(resample_dir, params.caas_full_perms == null ? 10 : (params.caas_full_perms as int), caas_config)  // 0 is a valid N: b_0 only

    emit:
        resample_subset = subset.subset
        fop_pairs = subset.fop_pairs.ifEmpty(file('NO_FOP_PAIRS'))
}
