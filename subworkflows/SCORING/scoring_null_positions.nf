#!/usr/bin/env nextflow
// scoring_null_positions.nf — Counters of the position-level permulation null for the position p-values.
// PhyloPhere | subworkflows/SCORING/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_NULL_POSITIONS: the integers that p.emp and p.emp_fact are made of, summed over the blocks of cycles of the
 *  null (null_position_stats.py documents what each counter is). The null is never held as a table of every
 *  (position, cycle); each block's table of scored positions is read once per stage.
 *
 *    1. CAAS_NULL_UNIVERSE   the positions the null scores (nd summed over the blocks) -> the universe U (tested set,
 *                            cleaned-background genes, observed positions) -> the propensity classes of p.emp_fact
 *    2. CAAS_NULL_COUNT      per block: k_emp and the per-class counts of each observed position
 *    3. CAAS_NULL_COUNTERS   the sums over the blocks: one row per observed position
 *
 *  Consumes:  null_positions.parquet and counts.json (CAAS_CORE_MERGE_REDUCE), the blocks of CAAS_CORE_MERGE, the roster,
 *             the observed filtered_discovery.tsv, the tested positions (background.output) and the cleaned background
 *             (each NO_* sentinel when absent)
 *  Produces:  null_counters.tsv, null_fit.tsv, null_universe_base.tsv (read by scoring_compute.R)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

process CAAS_NULL_UNIVERSE {
    tag "caas_null_universe"
    label 'process_medium'

    input:
    path null_positions
    path postproc_file
    path tested_positions   // caastools background.output or a NO_* sentinel
    path universe           // cleaned_background_main.txt or a NO_* sentinel
    path roster
    path counts
    path blocks, stageAs: 'blk_*'

    output:
    path "universe/null_classes.parquet", emit: classes
    path "universe/null_fit.tsv",         emit: fit
    path "universe/null_universe_base.tsv", emit: base

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local/src"
    def core_dir   = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    def tested_arg   = tested_positions.name.startsWith('NO_') ? '' : "--tested ${tested_positions}"
    def universe_arg = universe.name.startsWith('NO_') ? '' : "--universe ${universe}"
    """
    ${ (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript' } ${local_dir}/scoring_caas_perms.R --blocks blk_*/matrices.rds --cycles ${roster} --validate-only
    cp ${local_dir}/null_position_stats.py .
    mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py
    ${py} null_position_stats.py universe \\
        --null-positions ${null_positions} --postproc ${postproc_file} --roster ${roster} --counts-json ${counts} \\
        ${tested_arg} ${universe_arg} --output-dir universe
    """
}

process CAAS_NULL_COUNT {
    tag "caas_null_count|${block.name}"
    label 'process_medium'

    input:
    path block                 // a block directory of CAAS_CORE_MERGE (positions.parquet)
    path classes
    path fit
    path postproc_file

    output:
    path "counts.parquet", emit: counts

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local/src"
    def core_dir   = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    cp ${local_dir}/null_position_stats.py .
    mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py
    ${py} null_position_stats.py count \\
        --positions ${block}/positions.parquet --classes ${classes} --fit ${fit} --postproc ${postproc_file} \\
        --output counts.parquet
    """
}

process CAAS_NULL_COUNTERS {
    tag "caas_null_counters"
    label 'process_low'

    input:
    path counts, stageAs: 'counts_*.parquet'   // the counts of every block
    path classes
    path fit
    path postproc_file

    output:
    path "null_counters.tsv", emit: counters

    script:
    def local_dir  = "${baseDir}/subworkflows/SCORING/local/src"
    def core_dir   = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core"
    def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    cp ${local_dir}/null_position_stats.py .
    mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py
    ${py} null_position_stats.py reduce \\
        --counts counts_*.parquet --classes ${classes} --fit ${fit} --postproc ${postproc_file} \\
        --output null_counters.tsv
    """
}

workflow SCORING_NULL_POSITIONS {
    take:
        null_positions   // value channel: null_positions.parquet
        null_blocks      // value channel: the list of the blocks' directories
        null_counts      // value channel: counts.json
        postproc         // value channel: filtered_discovery.tsv
        tested_positions // value channel: background.output or a NO_* sentinel
        universe         // value channel: cleaned_background_main.txt or a NO_* sentinel
        roster           // value channel: cycle_roster.txt

    main:
        def u = CAAS_NULL_UNIVERSE(null_positions, postproc, tested_positions, universe, roster, null_counts, null_blocks)
        // The single-item outputs are queue channels: made value channels so that every block's task can read them.
        def classes_bc = u.classes.collect().map { items -> items[0] }
        def fit_bc     = u.fit.collect().map { items -> items[0] }
        def postproc_bc = postproc.collect().map { items -> items[0] }
        def per_block = CAAS_NULL_COUNT(null_blocks.flatMap { blocks -> blocks }, classes_bc, fit_bc, postproc_bc)
        def c = CAAS_NULL_COUNTERS(per_block.counts.collect(), classes_bc, fit_bc, postproc_bc)

    emit:
        counters = c.counters
        fit      = u.fit
        base     = u.base
}
