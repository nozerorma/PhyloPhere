#!/usr/bin/env nextflow

/*
 * CAAS permulation-excess null
 * ────────────────────────────
 * Builds a genome-wide *excess* null for CAAS FCS pathway enrichment:
 *   1. SUBSET_RESAMPLE_PERMS  — take the first N (caas_full_perms) permuted
 *      labelings from the resample output (drop b_0, the real labeling).
 *   2. BOOTSTRAP_PERMS        — full-pool bootstrap (no --discovery) with
 *      export_perm_discovery ON → per-gene per-cycle discovery rows.
 *   3. CONCAT_PERM_DISCOVERY  — stitch into one perm_discovery.tab.
 *   4. CAAS_PERMS_DISAMBIGUATE — load ASR once / replay N labelings
 *      (disambiguation_perms_main.py) → caas_perm_scores.tsv.
 *   5. CAAS_PERMS_AGGREGATE   — genes×N null matrices (scoring_caas_perms.R)
 *      → caas_perms.rds (corStat_byrank: global_asr/top_asr/bottom_asr).
 *
 * The aggregate RDS feeds the existing FCS p.perm path (fcs_enrich.R), giving
 * the CAAS scoring FCS report a permulation-corrected p.perm — exactly like RER.
 * See docs/CAAS_PERMULATION_EXCESS.md.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

// ── 1. Subset the resample to N permuted labelings (drop b_0) ────────────────
process SUBSET_RESAMPLE_PERMS {
    tag "caas_perms_subset|N=${n_perms}"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'resample_perms.tab'

    input:
    path resample_dir
    val  n_perms
    val  seed

    output:
    path "resample_perms.tab", emit: subset
    path "fop_pairs.tsv", emit: fop_pairs, optional: true

    script:
    """
    # Concatenate all resample cycles (-L: the resample dir is staged as a symlink,
    # so follow it). Drop the original labeling (b_0), then take a SEEDED RANDOM
    # sample of N labelings (not the first N — those are not exchangeable with an
    # arbitrary null draw). gawk's srand(seed) is deterministic, so the same seed
    # reproduces the same subset. We tag each candidate with a seeded rand() key,
    # sort by it, write to a file, then awk-limit by reading that FILE and exit-ing
    # early — avoids the SIGPIPE that `sort | head` triggers under pipefail.
    FOP_TAB=""
    if [ -d "${resample_dir}" ]; then
        FOP_TAB=\$(find -L ${resample_dir} -name 'fop_labelings.tab' | head -n 1)
    fi

    if [ -n "\$FOP_TAB" ]; then
        # ── FOP mirror: labelings are "<base>~H<m>". Sample N distinct BASE
        #    cycles and keep ALL their hypothesis rows together, plus the
        #    matching fop_pairs.tsv rows (PSS weights for domain pooling).
        FOP_PAIRS=\$(find -L ${resample_dir} -name 'fop_pairs.tsv' | head -n 1)
        awk -F'\\t' 'NF>=3 && \$1!="b_0"' "\$FOP_TAB" > candidates.tab
        awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' candidates.tab | sort -u > base_all.txt
        awk -v seed=${seed} 'BEGIN{srand(seed)} {print rand()"\\t"\$0}' base_all.txt \\
            | sort -k1,1g \\
            | awk -v n=${n_perms} '{sub(/^[^\\t]*\\t/,""); print; if(++c>=n) exit}' > keep_base.txt
        awk -F'\\t' 'NR==FNR{k[\$1]=1; next} {b=\$1; sub(/~.*/,"",b); if(b in k) print}' \\
            keep_base.txt candidates.tab > resample_perms.tab
        if [ -n "\$FOP_PAIRS" ]; then
            awk -F'\\t' 'NR==FNR{k[\$1]=1; next} FNR==1{if(!seen){print; seen=1}; next} (\$1 in k){print}' \\
                keep_base.txt "\$FOP_PAIRS" > fop_pairs.tsv
        fi
        n=\$(awk -F'\\t' '{b=\$1; sub(/~.*/,"",b); print b}' resample_perms.tab | sort -u | wc -l)
        rows=\$(wc -l < resample_perms.tab)
        avail=\$(wc -l < base_all.txt)
        echo "[caas_perms] FOP mirror seed=${seed}: \$n of \$avail base cycles (\$rows hypothesis labelings)"
        if [ "\$rows" -eq 0 ]; then echo "ERROR: no FOP labelings selected" >&2; exit 1; fi
        exit 0
    fi

    if [ -d "${resample_dir}" ]; then
        cat \$(find -L ${resample_dir} -name 'resample_*.tab' | sort) > all_resamples.tab
    elif [ -f "${resample_dir}" ]; then
        cat ${resample_dir} > all_resamples.tab
    else
        echo "ERROR: resample input not found: ${resample_dir}" >&2; exit 1
    fi
    awk -F'\\t' 'NF>=3 && \$1!="b_0"' all_resamples.tab > candidates.tab
    awk -v seed=${seed} 'BEGIN{srand(seed)} {print rand()"\\t"\$0}' candidates.tab \\
        | sort -k1,1g > shuffled.tab
    awk -v n=${n_perms} '{sub(/^[^\\t]*\\t/,""); print; if(++c>=n) exit}' \\
        shuffled.tab > resample_perms.tab
    n=\$(wc -l < resample_perms.tab)
    avail=\$(wc -l < candidates.tab)
    echo "[caas_perms] seed=${seed}: sampled \$n of \$avail permuted labelings (requested ${n_perms})"
    if [ "\$n" -eq 0 ]; then echo "ERROR: no permuted labelings selected" >&2; exit 1; fi
    """
}

// ── 2. Full-pool bootstrap with perm-discovery export (no --discovery) ───────
process BOOTSTRAP_PERMS {
    tag "$alignmentID"
    label 'process_boot'
    publishDir path: "${params.outdir}/caas_permulation/perm_disc", mode: 'copy', overwrite: true, pattern: '*.bootstrap.discovery.output'

    input:
    tuple val(alignmentID), path(alignmentFile), path(resampledPath)
    file caas_config

    output:
    tuple val(alignmentID), file("${alignmentID}.bootstrap.discovery.output"), emit: perm_discovery, optional: true

    script:
    def pairArgs = """
# multi_hypothesis mode passes a directory of traitfile_H*.tab (all share the
# same pair count K) — resolve it to one .tab so the awk reads a real file.
if [ -d "${caas_config}" ]; then
    _cfg_file=\$(find -L ${caas_config} -type f -name '*.tab' | head -n 1)
else
    _cfg_file="${caas_config}"
fi
n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$_cfg_file" | sort -nu | wc -l | tr -d ' ')
_max_conserved=\$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}')
_max_bg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_bg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
"""
    def ct_bin = (params.use_singularity || params.use_apptainer) ? "/usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct" : "$baseDir/subworkflows/CT/local/ct"
    """
    # Full-pool perms bootstrap now runs the vectorized CAAS kernel
    # (modules/boot_vec.py): the per-cycle CAAS test is a Level-3 BLAS matmul and
    # the perm_discovery rows are materialized only for the sparse hits. Size
    # OpenBLAS to the task allocation so those matmuls parallelize; pin MKL/NUMEXPR
    # (unused) and OpenMP to 1. Per-stage only — RER/R stages keep multithreaded
    # BLAS, so this must never go in nextflow.config env{}. See docs/CAAS_PERMULATION_RUNTIME.md.
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    ${pairArgs}
    # No --discovery → full position pool; export_perm_discovery forced ON.
    ${ct_bin} bootstrap \\
        -a ${alignmentFile} \\
        -t ${caas_config} \\
        -s ${resampledPath} \\
        -o ${alignmentID}.bootstraped.output \\
        --fmt ${params.ali_format} \\
        --patterns ${params.patterns} \\
        ${params.miss_pair ? '--miss_pair' : ''} \\
        ${params.caap_mode ? '--caap_mode' : ''} \\
        --export_perm_discovery ${alignmentID}.bootstrap.discovery.output \\
        --max_conserved \$_max_conserved \\
        --max_bg_gaps \$_max_bg_gaps \\
        --max_fg_gaps \$_max_fg_gaps \\
        --max_gaps \$_max_gaps \\
        --max_bg_miss \$_max_bg_miss \\
        --max_fg_miss \$_max_fg_miss \\
        --max_miss \$_max_miss
    """
}

// ── 2b. Batched full-pool bootstrap with perm-discovery export ───────────────
process BOOTSTRAP_PERMS_BATCHED {
    tag "$batchID (${batchSize} genes)"
    label 'process_boot_batched'
    publishDir path: "${params.outdir}/caas_permulation/perm_disc", mode: 'copy', overwrite: true, pattern: '*.bootstrap.discovery.output'

    input:
    tuple val(batchID), val(batchSize), val(batchManifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(resampledPath)
    file caas_config

    output:
    path("*.bootstraped.output"), emit: bootstrap_out, optional: true
    path("*.bootstrap.discovery.output"), emit: perm_discovery, optional: true

    script:
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct"
        : "$baseDir/subworkflows/CT/local/ct"
    """
    # Pin BLAS/OpenMP to 1 thread: the concurrency comes from the worker pool
    # running multiple genes in parallel; multi-threading would oversubscribe.
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

    # multi_hypothesis mode passes a directory of traitfile_H*.tab (same K) —
    # resolve to one .tab so awk reads a real file, not a directory.
    if [ -d "${caas_config}" ]; then
        _cfg_file=\$(find -L ${caas_config} -type f -name '*.tab' | head -n 1)
    else
        _cfg_file="${caas_config}"
    fi
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$_cfg_file" | sort -nu | wc -l | tr -d ' ')
    _max_conserved=\$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}')
    _max_bg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
    _max_fg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
    _max_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
    _max_bg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
    _max_fg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
    _max_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')

    declare -a extra_opts=()
    extra_opts+=(--patterns "${params.patterns}")
    if [ "${params.miss_pair}" = "true" ]; then extra_opts+=(--miss_pair); fi
    if [ "${params.caap_mode}" = "true" ]; then extra_opts+=(--caap_mode); fi
    extra_opts+=(--max_conserved \$_max_conserved)
    extra_opts+=(--max_bg_gaps \$_max_bg_gaps)
    extra_opts+=(--max_fg_gaps \$_max_fg_gaps)
    extra_opts+=(--max_gaps \$_max_gaps)
    extra_opts+=(--max_bg_miss \$_max_bg_miss)
    extra_opts+=(--max_fg_miss \$_max_fg_miss)
    extra_opts+=(--max_miss \$_max_miss)

    echo "\${extra_opts[@]}" > .ct_bootstrap_batch_args

    bash $baseDir/subworkflows/CT/local/scripts/run_ct_bootstrap_batch.sh \\
        --batch-id ${batchID} \\
        --manifest ${batchID}.manifest.tsv \\
        --caas-config ${caas_config} \\
        --resampled-path ${resampledPath} \\
        --workers ${task.cpus} \\
        --ali-format ${params.ali_format} \\
        --ct-bin ${ctBinary} \\
        --progress-log 0 \\
        --export-groups 0 \\
        --export-perm-discovery 1 \\
        --extra-args-file .ct_bootstrap_batch_args
    """
}

// ── 4. Null-mode disambiguation: load ASR once, replay N labelings ───────────
process CAAS_PERMS_DISAMBIGUATE {
    tag "caas_perms_disambiguate"
    label 'process_resample'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'gene_cycle_scores.tsv'

    input:
    path "perm_disc/*"
    path resample_subset
    path tree_file
    path fop_pairs   // fop_pairs.tsv (FOP mirror) or NO_FOP_PAIRS sentinel
    path gene_lengths // gene_ensembl_file (Gap B CT_POSTPROC filter) or NO_FILE

    output:
    path "gene_cycle_scores.tsv",    emit: gene_cycle_scores
    path "perm_pos_pval.tsv",        emit: pos_pval
    path "perm_pos_sample.tsv",      emit: pos_sample
    path "perm_pos_quantiles.tsv",   emit: pos_quantiles
    path "perm_pos_detail",          emit: pos_detail   // dir: one gz shard per gene

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def align_dir = params.alignment
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def workers = task.cpus ?: 1
    def max_tasks_per_child = params.ct_disambig_max_tasks_per_child ?: 50
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    // Gap B: mirror the observed CT_POSTPROC cluster + gene filters on the null
    // per-cycle CAAS pool. Opt-in via params.caas_perms_postproc (default true);
    // needs the gene annotation file for the extreme-gene density test.
    def _pp_raw = params.containsKey('caas_perms_postproc') ? params.caas_perms_postproc : true
    def _pp_on = (_pp_raw instanceof Boolean) ? _pp_raw : !(_pp_raw?.toString()?.toLowerCase() in ['false', '0', 'no'])
    def do_postproc = _pp_on && !(gene_lengths.name =~ /^NO_/)
    def postproc_args = do_postproc ? "--postproc-filter --gene-lengths ${gene_lengths} --clust-minlen ${params.filter_minlen} --clust-maxcaas ${params.filter_maxcaas} --gene-filter-mode ${params.gene_filter_mode} --iqr-multiplier ${params.iqr_multiplier} --extreme-percentile ${params.extreme_threshold}" : ""
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    mkdir -p caas_perms_out
    ${run} ./disambiguation_perms_main.py \\
        --alignment-dir ${align_dir} \\
        --tree ${tree_file} \\
        --perm-discovery perm_disc \\
        --resample-dir . \\
        --output-dir caas_perms_out \\
        ${fop_pairs.name =~ /^NO_/ ? '' : "--fop-pairs ${fop_pairs}"} ${postproc_args} \\
        --asr-model ${params.ct_disambig_asr_model} \\
        --convergence-mode ${params.ct_disambig_convergence_mode} \\
        --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
        --workers ${workers} \\
        --max-tasks-per-child ${max_tasks_per_child} \\
        --asr-cache-dir ${asr_cache_dir} \\
        ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
        ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    cp caas_perms_out/gene_cycle_scores.tsv gene_cycle_scores.tsv
    cp caas_perms_out/perm_pos_pval.tsv perm_pos_pval.tsv
    cp caas_perms_out/perm_pos_sample.tsv perm_pos_sample.tsv
    cp caas_perms_out/perm_pos_quantiles.tsv perm_pos_quantiles.tsv
    cp -R caas_perms_out/perm_pos_detail perm_pos_detail
    """
}

// ── 5. Aggregate to genes×N null matrices → caas_perms.rds ────────────────────
process CAAS_PERMS_AGGREGATE {
    tag "caas_perms_aggregate"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'caas_perms.rds'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'perm_pos_*.tsv'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'perm_pos_detail'

    input:
    path gene_cycle_scores
    path perm_pos_pval, stageAs: 'input_perm_pos_pval.tsv'
    path perm_pos_sample, stageAs: 'input_perm_pos_sample.tsv'
    path perm_pos_quantiles, stageAs: 'input_perm_pos_quantiles.tsv'
    path perm_pos_detail, stageAs: 'input_perm_pos_detail'   // dir: one gz shard per gene
    path universe

    output:
    path "caas_perms.rds",         emit: perms
    path "perm_pos_pval.tsv",      emit: pos_pval
    path "perm_pos_sample.tsv",    emit: pos_sample
    path "perm_pos_quantiles.tsv", emit: pos_quantiles
    path "perm_pos_detail",        emit: pos_detail
    path "gene_cycle_scores.tsv",  emit: gene_cycle_scores   // pass-through of the staged input; already published from DISAMBIGUATE

    script:
    def local_dir = "${baseDir}/subworkflows/SCORING/local"
    def universe_arg = universe.name != 'NO_FILE' ? "--universe ${universe}" : ""
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    cp ${perm_pos_pval} perm_pos_pval.tsv
    cp ${perm_pos_sample} perm_pos_sample.tsv
    cp ${perm_pos_quantiles} perm_pos_quantiles.tsv
    cp -R ${perm_pos_detail} perm_pos_detail

    ${run} ${local_dir}/src/scoring_caas_perms.R \\
        --gene-cycle-scores ${gene_cycle_scores} \\
        ${universe_arg} \\
        --output caas_perms.rds
    """
}

// ── 6. Rebuild the null from an existing detail file (no ASR replay) ─────────
// The gene-level null must hold the SAME statistic as the observed gene score, or
// the FCS p.perm compares two different quantities. A caas_perms.rds imported from
// a previous run is a cached artifact with no such guarantee — it was built with
// whatever formula was current then.
//
// Re-deriving it does NOT need the expensive part: perm_pos_detail/ (one gz shard
// per gene) already holds every (Gene, cycle, Position, caap_group,
// asr_path_score, n_detected, ct, cb) row, so only pass B (the aggregation) has
// to re-run. That is minutes against the hours of CAAS_PERMS_DISAMBIGUATE's ASR
// replay. A legacy single perm_pos_detail.tsv.gz file also still works.
//
// Deliberately reuses gene_wrapper.py's own aggregation via
// reaggregate_perm_scores.py rather than reimplementing: the gene statistic is
// F(max)^n over heavily tied values, so a 1e-16 difference in how the per-position
// sum accumulates can flip a tie and the ^n amplifies it.
process CAAS_PERMS_REBUILD {
    tag "caas_perms_rebuild"
    label 'process_medium'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true,
               pattern: '{caas_perms.rds,gene_cycle_scores.tsv,perm_pos_sample.tsv,perm_pos_quantiles.tsv}'

    input:
    path perm_pos_detail, stageAs: 'input_perm_pos_detail'   // dir (current) or legacy .tsv.gz file
    path universe

    output:
    path "caas_perms.rds",         emit: perms
    path "gene_cycle_scores.tsv",  emit: gene_cycle_scores
    path "perm_pos_sample.tsv",    emit: pos_sample,    optional: true
    path "perm_pos_quantiles.tsv", emit: pos_quantiles, optional: true

    script:
    def disambig_local = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def scoring_local  = "${baseDir}/subworkflows/SCORING/local"
    // Sentinel convention is a NO_* basename (NO_FILE here, NO_BACKGROUND from
    // SCORING's resolved_background) — match the prefix, not one literal.
    def universe_arg   = universe.name.startsWith('NO_') ? "" : "--universe ${universe}"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript'  : 'Rscript'
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    # reaggregate_perm_scores.py imports src.utils.gene_wrapper relative to its own
    # location, so stage the disambiguation local tree exactly as CAAS_PERMS_DISAMBIGUATE does.
    cp -R ${disambig_local}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true

    ${py} ./reaggregate_perm_scores.py \\
        --detail input_perm_pos_detail \\
        --output-dir .

    ${rs} ${scoring_local}/src/scoring_caas_perms.R \\
        --gene-cycle-scores gene_cycle_scores.tsv \\
        ${universe_arg} \\
        --output caas_perms.rds
    """
}

// ── Subworkflow: prep (subset + full-pool export) — runs INSIDE ct.nf where the
//    sliced per-gene alignments (align_tuple) + resample_dir are available. ────
workflow CAAS_PERMS_PREP {
    take:
        align_tuple        // Channel<tuple(id, alignmentFile)>
        caas_config        // path (trait file)
        resample_dir       // path (resample_*.tab directory)

    main:
        def bootstrapBatchSize = (params.ct_bootstrap_batch_size ?: 1) as int
        def subset = SUBSET_RESAMPLE_PERMS(resample_dir, params.caas_full_perms ?: 10, params.seed ?: 1998)
        def boot_in = align_tuple
            .map { id, f -> tuple(id, f) }
            .combine(subset.subset)
            .map { id, f, sub -> tuple(id, f, sub) }

        def discovery_files
        if (bootstrapBatchSize > 1) {
            def bootstrapBatchCounter = 0
            def bootstrap_batches = boot_in
                .toSortedList({ a, b -> a[0] <=> b[0] })
                .flatMap()
                .collate(bootstrapBatchSize)
                .map { batch ->
                    def batchID = sprintf('bootstrap_perms_batch_%05d', ++bootstrapBatchCounter)
                    def manifestText = batch.collect { row -> "${row[0]}\t${row[1].name}\tNO_FILE" }.join('\n') + '\n'
                    def alignmentFiles = batch.collect { row -> row[1] }.unique { file -> file.name }
                    def resampled = batch[0][2]
                    tuple(batchID, batch.size(), manifestText, alignmentFiles, resampled)
                }
            def boot = BOOTSTRAP_PERMS_BATCHED(bootstrap_batches, caas_config)
            discovery_files = boot.perm_discovery.flatten().collect()
        } else {
            def boot = BOOTSTRAP_PERMS(boot_in, caas_config)
            discovery_files = boot.perm_discovery.map { id, f -> f }.collect()
        }
    emit:
        perm_discovery = discovery_files
        resample_subset = subset.subset
        fop_pairs = subset.fop_pairs.ifEmpty(file('NO_FOP_PAIRS'))
}

// ── Subworkflow: build the null matrices — runs in main.nf after CT ───────────
workflow CAAS_PERMULATION {
    take:
        perm_discovery     // path perm_discovery.tab
        resample_subset    // path resample_perms.tab
        tree_file          // path species tree
        universe           // path cleaned_background or NO_FILE
        fop_pairs          // path fop_pairs.tsv (FOP mirror) or NO_FOP_PAIRS
        gene_lengths       // path gene_ensembl_file (Gap B CT_POSTPROC filter) or NO_FILE
        asr_ready          // gate: emits once the live CT_DISAMBIGUATION has
                           // finished writing the shared ASR cache. A 'NO_GATE'
                           // sentinel when ASR is precomputed / disambiguation
                           // is off (no wait needed then).

    main:
        // CAAS_PERMS_DISAMBIGUATE ("load precomputed ASR once, replay N
        // labelings") has NO compute-on-miss — it only reads --asr-cache-dir. In
        // asr_mode=compute the live CT_DISAMBIGUATION_RUN populates that dir
        // lazily and there is no other dependency edge between the two, so on a
        // fast run the replay reads the cache before it is written and the whole
        // permulation null comes out empty, silently. Hold the tree input until
        // asr_ready emits so the replay always runs after the cache is complete.
        def gated_tree = tree_file
            .combine(asr_ready)
            .map { t, _ready -> t }
        def scores = CAAS_PERMS_DISAMBIGUATE(perm_discovery, resample_subset, gated_tree, fop_pairs, gene_lengths)
        def agg = CAAS_PERMS_AGGREGATE(scores.gene_cycle_scores, scores.pos_pval, scores.pos_sample,
                                       scores.pos_quantiles, scores.pos_detail, universe)

    emit:
        perms              = agg.perms
        pos_pval           = agg.pos_pval        // LOO null_pvalue_boot per (gene,position,scheme)
        pos_sample         = agg.pos_sample      // cycle-stratified sample for distribution plots
        pos_quantiles      = agg.pos_quantiles   // per (cycle,scheme) distribution shape
        pos_detail         = agg.pos_detail      // full per-cycle detail (sharded dir); re-scoring needs no ASR replay
        gene_cycle_scores  = agg.gene_cycle_scores // genes x cycles raw scores; feeds the report's FPR calibration figure
}
