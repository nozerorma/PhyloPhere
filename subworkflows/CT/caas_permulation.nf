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

    input:
    path resample_dir
    val  n_perms

    output:
    path "resample_perms.tab", emit: subset

    script:
    """
    # Concatenate all resample cycles (-L: the resample dir is staged as a symlink,
    # so follow it). Drop the original labeling (b_0) and keep the first N. awk does
    # the limiting by reading a FILE and exit-ing early — avoids the SIGPIPE that
    # `... | head` triggers under pipefail.
    cat \$(find -L ${resample_dir} -name 'resample_*.tab' | sort) > all_resamples.tab
    awk -v n=${n_perms} -F'\\t' 'NF>=3 && \$1!="b_0"{print; if(++c>=n) exit}' \\
        all_resamples.tab > resample_perms.tab
    n=\$(wc -l < resample_perms.tab)
    echo "[caas_perms] selected \$n permuted labelings (requested ${n_perms})"
    if [ "\$n" -eq 0 ]; then echo "ERROR: no permuted labelings selected" >&2; exit 1; fi
    """
}

// ── 2. Full-pool bootstrap with perm-discovery export (no --discovery) ───────
process BOOTSTRAP_PERMS {
    tag "$alignmentID"
    label 'process_boot'

    input:
    tuple val(alignmentID), path(alignmentFile), path(resampledPath)
    file caas_config

    output:
    tuple val(alignmentID), file("${alignmentID}.bootstrap.discovery.output"), emit: perm_discovery, optional: true

    script:
    def pairArgs = """
n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' ${caas_config} | sort -nu | wc -l | tr -d ' ')
_max_conserved=\$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}')
_max_bg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_bg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
"""
    def ct_bin = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh ct' : "$baseDir/subworkflows/CT/local/ct"
    """
    # Pin BLAS/OpenMP to 1 thread: the CAAS kernel is scalar scipy.stats (no gemm),
    # so library auto-threading is pure overhead and oversubscribes under maxForks.
    # (Not global — R/RER stages keep their multithreaded BLAS.) See docs/CAAS_PERMULATION_RUNTIME.md.
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
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

// ── 3. Concatenate per-gene perm-discovery into one table (single header) ─────
process CONCAT_PERM_DISCOVERY {
    tag "caas_perms_concat"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'perm_discovery.tab'

    input:
    path discovery_files

    output:
    path "perm_discovery.tab", emit: perm_discovery

    script:
    """
    set -euo pipefail
    files=( ${discovery_files} )
    first=\${files[0]}
    head -n 1 "\$first" > perm_discovery.tab
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> perm_discovery.tab
    done
    echo "[caas_perms] concatenated \$(( \$(wc -l < perm_discovery.tab) - 1 )) discovery rows"
    """
}

// ── 4. Null-mode disambiguation: load ASR once, replay N labelings ───────────
process CAAS_PERMS_DISAMBIGUATE {
    tag "caas_perms_disambiguate"
    label 'process_resample'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'caas_perm_scores.tsv'

    input:
    path perm_discovery
    path resample_subset
    path tree_file

    output:
    path "caas_perm_scores.tsv", emit: perm_scores

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def align_dir = params.alignment
    def asr_cache_dir = params.ct_disambig_asr_cache_dir ?: ''
    def taxid_mapping = params.tax_id ?: ''
    def ensembl_file = params.gene_ensembl_file ?: ''
    def workers = params.ct_disambig_workers ?: (task.cpus ?: 1)
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    # Pin BLAS/OpenMP to 1 thread BEFORE the mp.Pool launches — each replay worker
    # inherits this, so the pool's parallelism is the workers, not workers×nproc BLAS
    # threads. The replay (path_scores) is pure-Python; nothing here benefits from BLAS.
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    cp -R ${local_dir}/* .
    find . -name '*.pyc' -delete 2>/dev/null || true
    mkdir -p caas_perms_out
    ${run} ./disambiguation_perms_main.py \\
        --alignment-dir ${align_dir} \\
        --tree ${tree_file} \\
        --perm-discovery ${perm_discovery} \\
        --resample-dir . \\
        --output-dir caas_perms_out \\
        --asr-model ${params.ct_disambig_asr_model} \\
        --convergence-mode ${params.ct_disambig_convergence_mode} \\
        --posterior-threshold ${params.ct_disambig_posterior_threshold} \\
        --workers ${workers} \\
        --asr-cache-dir ${asr_cache_dir} \\
        ${taxid_mapping ? "--taxid-mapping ${taxid_mapping}" : ''} \\
        ${ensembl_file ? "--ensembl-genes-file ${ensembl_file}" : ''}
    cp caas_perms_out/caas_perm_scores.tsv caas_perm_scores.tsv
    """
}

// ── 5. Aggregate to genes×N null matrices → caas_perms.rds ────────────────────
process CAAS_PERMS_AGGREGATE {
    tag "caas_perms_aggregate"
    label 'process_low'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'caas_perms.rds'
    publishDir path: "${params.outdir}/caas_permulation", mode: 'copy', overwrite: true, pattern: 'perm_pos_*.tsv'

    input:
    path perm_scores
    path universe

    output:
    path "caas_perms.rds",      emit: perms
    path "perm_pos_pval.tsv",   emit: pos_pval     // recovery p-value per (gene,position,scheme)
    path "perm_pos_sample.tsv", emit: pos_sample   // capped per-scheme sample for distribution plots

    script:
    def local_dir = "${baseDir}/subworkflows/SCORING/local"
    def universe_arg = universe.name != 'NO_FILE' ? "--universe ${universe}" : ""
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'
    """
    ${run} ${local_dir}/scoring_caas_perms.R \\
        --perm-scores ${perm_scores} \\
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
        def subset = SUBSET_RESAMPLE_PERMS(resample_dir, params.caas_full_perms ?: 10)
        def boot_in = align_tuple
            .map { id, f -> tuple(id, f) }
            .combine(subset.subset)
            .map { id, f, sub -> tuple(id, f, sub) }
        def boot = BOOTSTRAP_PERMS(boot_in, caas_config)
        def discovery_files = boot.perm_discovery.map { id, f -> f }.collect()
        def concat = CONCAT_PERM_DISCOVERY(discovery_files)

    emit:
        perm_discovery = concat.perm_discovery
        resample_subset = subset.subset
}

// ── Subworkflow: build the null matrices — runs in main.nf after CT ───────────
workflow CAAS_PERMULATION {
    take:
        perm_discovery     // path perm_discovery.tab
        resample_subset    // path resample_perms.tab
        tree_file          // path species tree
        universe           // path cleaned_background or NO_FILE

    main:
        def scores = CAAS_PERMS_DISAMBIGUATE(perm_discovery, resample_subset, tree_file)
        def agg = CAAS_PERMS_AGGREGATE(scores.perm_scores, universe)

    emit:
        perms       = agg.perms
        perm_scores = scores.perm_scores   // per-(gene,cycle,position,scheme) null asr_path_score
        pos_pval    = agg.pos_pval          // lean recovery p-value per (gene,position,scheme)
        pos_sample  = agg.pos_sample        // lean capped sample for distribution plots
}
