#!/usr/bin/env nextflow
// scoring_compute.nf — Position-level and gene-level CAAS scores, joined with FADE, RER and accumulation evidence.
// PhyloPhere | subworkflows/SCORING/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_COMPUTE: scores the observed positions and genes of the filtered discovery
 *  and integrates the other tools' evidence per gene.
 *
 *  observed_core_scores.py first computes the position and gene CAAS scores with the
 *  code that also scores the permulation null (src/core/scores.py, copied from
 *  CT_DISAMBIGUATION); scoring_compute.R then adds the position and gene permulation
 *  p-values, joins FADE, RERConverge and accumulation, and writes the tables, the
 *  ranked slices and the enrichment curves. It runs once on the full postproc pool;
 *  direction is carried by the side column.
 *
 *  Consumes:  postproc_file (filtered_discovery.tsv, mandatory); optional FADE gene
 *             summaries and per-site tables, RERConverge summary, accumulation results,
 *             hypothesis pairs, CAAS permulation null (caas_perms.rds and
 *             perm_pos_cycle_caas.tsv.gz). Each optional input is a staged
 *             NO_* sentinel file when absent.
 *  Produces:  position_scores.tsv, gene_scores.tsv, gene_correlations.tsv,
 *             fcs_stats.tsv (plus fcs_stats_{rer,fade,accum}.tsv when that evidence
 *             exists), gene_lists/ and position_lists/ (12 slice TSVs each: top,
 *             bottom, global × 25, 10, 5, 1%; position slices keyed Gene, Position,
 *             CAAS_score), gene_threshold_enrichment.tsv, pos_threshold_enrichment.tsv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Scoring ────────────────────────────────────────────────────────────────────

process SCORING_COMPUTE {
    tag "scoring_compute|${params.traitname ?: 'unknown_trait'}"
    label 'error_retry'

    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: 'gene_lists'
    publishDir path: "${params.outdir}/scoring",
               mode: 'copy', overwrite: true,
               pattern: 'position_lists'

    input:
    path postproc_file
    path fade_summary_top
    path fade_summary_bottom
    path fade_site_top
    path fade_site_bot
    path rer_summary
    path accum_files
    path hypotheses_pairs   // optional: contrast_hypotheses_pairs.tsv (FOP); NO_HYP_PAIRS sentinel otherwise
    path caas_perms         // optional: caas_perms.rds, the gene-level permulation null; NO_FILE sentinel otherwise
    path caas_pos_cycle_caas // optional: perm_pos_cycle_caas.tsv.gz, the position-level null for p.emp; NO_FILE sentinel otherwise

    output:
    path "position_scores.tsv",                              emit: position_scores
    path "gene_scores.tsv",                                  emit: gene_scores
    path "fcs_stats.tsv",                                    emit: fcs_stats
    path "fcs_stats_rer.tsv",                          optional: true, emit: fcs_stats_rer
    path "fcs_stats_fade.tsv",                         optional: true, emit: fcs_stats_fade
    path "fcs_stats_accum.tsv",                        optional: true, emit: fcs_stats_accum
    path "gene_correlations.tsv",                            emit: gene_correlations
    path "gene_lists",                                 optional: true, emit: gene_lists
    path "position_lists",                             optional: true, emit: position_lists
    path "gene_threshold_enrichment.tsv",             optional: true, emit: gene_threshold_enrichment
    path "pos_threshold_enrichment.tsv",              optional: true, emit: pos_threshold_enrichment

    script:
    def local_dir       = "${baseDir}/subworkflows/SCORING/local/src"
    def core_dir        = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local/src/core"
    def top_pct         = params.scoring_position_top_pct    ?: 0.10
    def g_top_pct       = params.scoring_gene_top_pct        ?: 0.10
    def accum_arg         = (accum_files instanceof List
                              ? (accum_files.size() == 1 && accum_files[0].name.startsWith('NO_') ? 'NO_ACCUM' : '.')
                              : (accum_files.name.startsWith('NO_') ? 'NO_ACCUM' : '.'))
    def fs_top_arg        = fade_site_top.name =~ /^NO_FADE_SITE_TOP/ ? 'NO_FADE_SITE_TOP' : "${fade_site_top}"
    def fs_bot_arg        = fade_site_bot.name =~ /^NO_FADE_SITE_BOT/ ? 'NO_FADE_SITE_BOT' : "${fade_site_bot}"
    def hp_arg            = hypotheses_pairs.name =~ /^NO_/ ? 'NO_HYP_PAIRS' : "${hypotheses_pairs}"
    def cp_arg            = caas_perms.name =~ /^NO_/ ? 'NO_FILE' : "${caas_perms}"
    def cpcc_arg          = caas_pos_cycle_caas.name =~ /^NO_/ ? 'NO_FILE' : "${caas_pos_cycle_caas}"
    def gene_perm_pooled  = params.scoring_gene_perm_pooled ?: false

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/scoring_compute.R ${local_dir}/aa_grouping.R ${local_dir}/observed_core_scores.py .
        mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py

        /usr/local/bin/_entrypoint.sh python3 observed_core_scores.py \
            --input '${postproc_file}' \
            --positions-out core_positions.tsv \
            --genes-out core_genes.tsv

        /usr/local/bin/_entrypoint.sh Rscript scoring_compute.R \
            --postproc       '${postproc_file}' \
            --core_positions core_positions.tsv \
            --core_genes     core_genes.tsv \
            --fade_top       '${fade_summary_top}' \
            --fade_bottom    '${fade_summary_bottom}' \
            --fade_site_top  '${fs_top_arg}' \
            --fade_site_bot  '${fs_bot_arg}' \
            --rer            '${rer_summary}' \
            --accum_dir      '${accum_arg}' \
            --hypotheses_pairs '${hp_arg}' \
            --caas_perms      '${cp_arg}' \
            --caas_pos_cycle_caas '${cpcc_arg}' \
            --gene_perm_pooled '${gene_perm_pooled}' \
            --p_emp_thr             ${params.scoring_p_emp_thr ?: 0.05} \
            --top_pct              ${top_pct} \
            --top25_pct            0.25 \
            --top5_pct             0.05 \
            --top1_pct             0.01 \
            --gene_top_pct         ${g_top_pct} \
            --gene_top25_pct       0.25 \
            --gene_top5_pct        0.05 \
            --gene_top1_pct        0.01
        """
    } else {
        """
        cp ${local_dir}/scoring_compute.R ${local_dir}/aa_grouping.R ${local_dir}/observed_core_scores.py .
        mkdir -p src/core && cp ${core_dir}/scores.py src/core/ && touch src/__init__.py src/core/__init__.py

        python3 observed_core_scores.py \
            --input '${postproc_file}' \
            --positions-out core_positions.tsv \
            --genes-out core_genes.tsv

        Rscript scoring_compute.R \
            --postproc       '${postproc_file}' \
            --core_positions core_positions.tsv \
            --core_genes     core_genes.tsv \
            --fade_top       '${fade_summary_top}' \
            --fade_bottom    '${fade_summary_bottom}' \
            --fade_site_top  '${fs_top_arg}' \
            --fade_site_bot  '${fs_bot_arg}' \
            --rer            '${rer_summary}' \
            --accum_dir      '${accum_arg}' \
            --hypotheses_pairs '${hp_arg}' \
            --caas_perms      '${cp_arg}' \
            --caas_pos_cycle_caas '${cpcc_arg}' \
            --gene_perm_pooled '${gene_perm_pooled}' \
            --p_emp_thr             ${params.scoring_p_emp_thr ?: 0.05} \
            --top_pct              ${top_pct} \
            --top25_pct            0.25 \
            --top5_pct             0.05 \
            --top1_pct             0.01 \
            --gene_top_pct         ${g_top_pct} \
            --gene_top25_pct       0.25 \
            --gene_top5_pct        0.05 \
            --gene_top1_pct        0.01
        """
    }
}
