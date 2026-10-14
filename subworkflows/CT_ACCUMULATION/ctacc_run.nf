#!/usr/bin/env nextflow
// ctacc_run.nf — Alignment variability, position aggregation and randomization for CAAS accumulation.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT_ACCUMULATION_AGGREGATE and CT_ACCUMULATION_RANDOMIZE: the two phases of
 *  local/main.py, and COMPUTE_ALIGNMENT_ENTROPY, which prepares the conservation input
 *  of the first one. Each script block has a container branch (entrypoint wrapper) and
 *  a plain branch that run the same command.
 *
 *  Consumes:  alignment directory, genomic-info TSV, traitfile, filtered_discovery.tsv
 *             and cleaned background list (CT_POSTPROC); for the randomization also
 *             the tested positions (background.output) and perm_pos_cycle_caas.tsv.gz and
 *             cycle_roster.txt of one CAAS_CORE_MERGE run (the null of the permulation type,
 *             the universe U of the draw pool of naive and cons_decile)
 *  Produces:  entropy_dir/ (<gene>.entropy.tsv), accumulation_global.csv,
 *             accumulation_<direction>_aggregated_results.csv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Alignment variability ────────────────────────────────────────────────────

// Valdar variability per column and gene (compute_alignment_entropy.py); feeds the conservation value of the aggregation
process COMPUTE_ALIGNMENT_ENTROPY {
    tag "auto-generate accumulation entropy"
    label 'process_medium'

    publishDir "${params.outdir}/core_inputs", mode: 'copy', overwrite: true

    input:
    path alignment_dir
    path taxid_tsv

    output:
    path "entropy_dir", emit: entropy_dir

    script:
    """
    python3 ${baseDir}/bin/compute_alignment_entropy.py \\
        --alignment-dir "${alignment_dir}" \\
        --output-dir entropy_dir \\
        --taxid-tsv "${taxid_tsv}"
    """
}

// ── Aggregation ──────────────────────────────────────────────────────────────

// One global position table (<prefix>_global.csv) for all directions
process CT_ACCUMULATION_AGGREGATE {
    tag "ct_accumulation_aggregate"
    label 'process_long_compute'

    publishDir path: "${params.outdir}/accumulation/aggregation", mode: 'copy', overwrite: true,
               pattern: '*.csv'

    input:
    val  alignment_dir
    path genomic_info
    path species_list
    path metadata_caas
    path bg_caas
    val  entropy_dir

    output:
    path "*_global.csv",    emit: global_csv

    script:
    def local_dir    = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def ali_fmt      = params.ali_format
    def out_pfx      = 'accumulation'
    def log_level    = 'INFO'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        /usr/local/bin/_entrypoint.sh python main.py \\
            --tool aggregate \\
            --alignment-dir "${alignment_dir}" \\
            --alignment-format '${ali_fmt}' \\
            --genomic-info '${genomic_info}' \\
            --species-list '${species_list}' \\
            --metadata-caas '${metadata_caas}' \\
            --bg-caas '${bg_caas}' \\
            --output-prefix '${out_pfx}' \\
            --entropy-dir '${entropy_dir}' \\
            --log-level '${log_level}'
        """
    } else {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        python main.py \\
            --tool aggregate \\
            --alignment-dir "${alignment_dir}" \\
            --alignment-format '${ali_fmt}' \\
            --genomic-info '${genomic_info}' \\
            --species-list '${species_list}' \\
            --metadata-caas '${metadata_caas}' \\
            --bg-caas '${bg_caas}' \\
            --output-prefix '${out_pfx}' \\
            --entropy-dir '${entropy_dir}' \\
            --log-level '${log_level}'
        """
    }
}

// ── Randomization ────────────────────────────────────────────────────────────

// One run per direction ('top', 'bottom', 'all'): the null and the per-gene empirical p-values
process CT_ACCUMULATION_RANDOMIZE {
    tag "ct_accumulation_randomize|${direction}"
    label 'process_long_compute'

    publishDir path: { "${params.outdir}/accumulation/${direction}/randomization" },
               mode: 'copy', overwrite: true,
               pattern: '*_aggregated_results.csv'

    input:
    val  direction
    path global_csv
    path caas_csv
    path background_positions   // caastools background.output (tested positions)
    path bg_caas_universe       // cleaned_background_main.txt (surviving genes)
    path null_cycle_caas        // perm_pos_cycle_caas.tsv.gz from CAAS_CORE_MERGE, the filtered null SCORING also uses (NO_ sentinel when absent):
                                 // the null of 'permulation', the universe U of the draw pool of 'naive' and 'cons_decile'
    path cycle_roster           // cycle_roster.txt of the same run (NO_ sentinel when absent): N of the permulation null,
                                 // counting the cycles that detected nothing (see randomize.py)

    output:
    val  direction,                           emit: direction
    path "*_aggregated_results.csv",          emit: results

    script:
    def local_dir    = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def out_pfx        = "accumulation_${direction}"
    // The eligible null pool is the tested positions of the surviving genes, restricted to the positions
    // the null scores plus the observed ones (universe U) when the null table is given. Every input
    // is optional at the channel level (a NO_ sentinel file when absent); randomize.py
    // falls back to the ungapped-column pool, with a warning, when the tested positions are missing,
    // and to the tested pool, logged, when the null table is missing.
    def bgpos_flag = (background_positions.name =~ /^NO_/) ? '' : "--background-positions '${background_positions}'"
    def bgcaas_flag = (bg_caas_universe.name =~ /^NO_/) ? '' : "--bg-caas '${bg_caas_universe}'"
    // The null table feeds every randomization type (the null of permulation, the universe U of naive and
    // cons_decile); the roster only matters to permulation, which takes N from it.
    def nullcycle_flag = (null_cycle_caas.name =~ /^NO_/) ? '' : "--null-cycle-caas '${null_cycle_caas}'"
    def roster_flag = (cycle_roster.name =~ /^NO_/) ? '' : "--cycle-roster '${cycle_roster}'"
    // 'all' maps to --change-side both, which keeps every row with a side other than 'none'
    def change_side_arg = (direction == 'all') ? 'both' : direction
    def rand_type    = params.accumulation_randomization_type ?: 'naive'
    def n_rands      = params.accumulation_n_randomizations   ?: 10000
    def log_level    = 'INFO'
    // Always pass task.cpus: without --workers randomize.py uses os.cpu_count(), which is
    // the core count of the node, not of the Slurm allocation.
    def workers_flag = "--workers ${task.cpus}"
    def seed_flag    = params.seed ? "--global-seed ${params.seed}" : '--global-seed 1998'

    if (params.use_singularity || params.use_apptainer) {
        """
        # The randomizations are split across a ProcessPoolExecutor of numpy workers
        # (random draws and bincount, no BLAS-heavy work). Math-library threads are pinned
        # to 1 per worker so that library threading does not multiply by the worker count.
        # Set in this process only, not in a global env scope.
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        /usr/local/bin/_entrypoint.sh python main.py \\
            --tool randomize \\
            --global-csv '${global_csv}' \\
            --caas-csv '${caas_csv}' \\
            --output-prefix '${out_pfx}' \\
            --randomization-type '${rand_type}' \\
            --n-randomizations ${n_rands} \\
            --change-side '${change_side_arg}' \\
            ${bgpos_flag} ${bgcaas_flag} \\
            ${nullcycle_flag} ${roster_flag} \\
            ${workers_flag} ${seed_flag} \\
            --log-level '${log_level}'
        """
    } else {
        """
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        find . -name '*.pyc' -delete 2>/dev/null || true

        python main.py \\
            --tool randomize \\
            --global-csv '${global_csv}' \\
            --caas-csv '${caas_csv}' \\
            --output-prefix '${out_pfx}' \\
            --randomization-type '${rand_type}' \\
            --n-randomizations ${n_rands} \\
            --change-side '${change_side_arg}' \\
            ${bgpos_flag} ${bgcaas_flag} \\
            ${nullcycle_flag} ${roster_flag} \\
            ${workers_flag} ${seed_flag} \\
            --log-level '${log_level}'
        """
    }
}
