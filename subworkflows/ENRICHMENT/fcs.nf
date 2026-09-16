#!/usr/bin/env nextflow

/*
 * FCS (Functional Class Scoring) report processes
 * ───────────────────────────────────────────────
 * Rank-based, threshold-free gene-set enrichment via the Wilcoxon-AUC test
 * (RERconverge::fastwilcoxGMTall) over the curated GMTs in subworkflows/ENRICHMENT/dat.
 * Each process renders 12.FCS_general_report.Rmd against a generic
 * stats TSV (gene + score_<ranking> + flag_<name> columns) and a universe file
 * (cleaned_background, no-signal genes floored to 0).
 *
 *   SCORING_FCS_REPORT : CAAS scoring (global/top/bottom + full cross-module flags)
 *   RER_FCS_REPORT     : RERconverge-specific report process (takes perms_file)
 *
 * FADE and Accumulation do not run their own FCS ranking: FADE's statistic is a
 * max over many sites and Accumulation has no permulation null, so neither
 * supports a reliable standalone significance test. They contribute as
 * cross-module corroboration flags on CAAS's/RER's leading edge instead, and
 * FADE additionally gets its own position-level group in posenrich.
 */

// ─────────────────────────────────────────────────────────────────────────────
// SCORING_FCS_REPORT
// ─────────────────────────────────────────────────────────────────────────────
process SCORING_FCS_REPORT {
    tag "scoring_fcs|${params.traitname ?: 'unknown_trait'}"
    label 'process_reporting'

    publishDir path: "${params.outdir}/fcs",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/fcs/fcs_results",
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**'

    input:
    path fcs_stats
    path universe
    path perms_file
    path gene_lists
    path enrich_file

    output:
    path "12.FCS_scoring_${params.traitname ?: 'unknown_trait'}.html", emit: report
    path "fcs_results/**",                       emit: fcs_results,      optional: true
    path "fcs_results/fcs_all_results.tsv",      emit: fcs_all_results,  optional: true
    path "fcs_results/fcs_leading_edge.tsv",     emit: fcs_leading_edge, optional: true
    path "fcs_results/fcs_leading_edge_composition.tsv", emit: fcs_leading_edge_composition, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def traitname = params.traitname ?: 'unknown_trait'
    def gmt_dir   = params.gmt_dir
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def fdr_permsum     = params.fcs_fdr_permsum     ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    // SCORING's own published gene_lists/ -- this IS the CAAS report, so
    // score_top/score_bottom here really are CAAS's gene_caas_score_top_all/
    // bottom_all (see 12.FCS_general_report.Rmd's gene_lists_dir param doc).
    def gene_lists_arg = (gene_lists.name =~ /^NO_/) ? 'NULL' : "'${gene_lists}'"
    def enrich_file_arg = (enrich_file.name =~ /^NO_/) ? 'NULL' : "'${enrich_file}'"
    def render = """
        rmarkdown::render(
            '12.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = 'Scoring_FCS_${traitname}',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                fdr_wilcoxon    = ${fdr_wilcoxon},
                fdr_lachenbruch = ${fdr_lachenbruch},
                fdr_permsum     = ${fdr_permsum},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${traitname}',
                perms_file    = '${perms_file}',
                gene_lists_dir = ${gene_lists_arg},
                enrich_file   = ${enrich_file_arg},
                seed          = '${params.seed ?: 1998}'
            ),
            output_file = '12.FCS_scoring_${traitname}.html'
        )
    """
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "${render}"
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "${render}"
        """
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RER_FCS_REPORT - RERconverge specific report process that takes perms_file
// ─────────────────────────────────────────────────────────────────────────────
process RER_FCS_REPORT {
    tag "rer_fcs|${report_label}"
    label 'process_reporting'

    publishDir path: { "${params.outdir}/${subpath.toLowerCase()}" },
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true, pattern: '*.html'
    publishDir path: { "${params.outdir}/${subpath.toLowerCase()}/fcs_results" },
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**'

    input:
    val  subpath
    path fcs_stats
    path universe
    val  report_label
    path perms_file
    path annot_file
    path enrich_file

    output:
    path "${report_label}.html",             emit: report
    path "fcs_results/**",                   emit: fcs_results,      optional: true
    path "fcs_results/fcs_all_results.tsv",  emit: fcs_all_results,  optional: true
    path "fcs_results/fcs_leading_edge.tsv", emit: fcs_leading_edge, optional: true
    path "fcs_results/fcs_leading_edge_composition.tsv", emit: fcs_leading_edge_composition, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def gmt_dir   = params.gmt_dir
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def fdr_permsum     = params.fcs_fdr_permsum     ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    def enrich_file_arg = (enrich_file.name =~ /^NO_/) ? 'NULL' : "'${enrich_file}'"
    def render = """
        rmarkdown::render(
            '12.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = '${report_label}',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                fdr_wilcoxon    = ${fdr_wilcoxon},
                fdr_lachenbruch = ${fdr_lachenbruch},
                fdr_permsum     = ${fdr_permsum},
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${params.traitname ?: "trait"}',
                perms_file    = '${perms_file}',
                annot_file    = '${annot_file}',
                enrich_file   = ${enrich_file_arg},
                seed          = '${params.seed ?: 1998}'
            ),
            output_file = '${report_label}.html'
        )
    """
    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        /usr/local/bin/_entrypoint.sh Rscript -e "${render}"
        """
    } else {
        """
        cp -R ${local_dir}/* .
        Rscript -e "${render}"
        """
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FCS_COMPUTE_BATCHED / FCS_CONCAT / FCS_COMPUTE
// ─────────────────────────────────────────────────────────────────────────────
// fcs_run_all() over every GMT database is the expensive step behind both
// report processes above (Wilcoxon-AUC + Lachenbruch + Path-Sum-Permulation for
// every score_<ranking> column, against every GMT). BH correction throughout
// fcs_enrich.R is scoped per database (fastwilcoxGMTall's own per-GMT BH, plus
// the explicit per-db p.adjust() calls for the one-sided recompute, Lachenbruch,
// and Path-Sum-Permulation - see that file's header), so splitting the GMT set
// across independent batched Nextflow tasks and row-concatenating their results
// is exact, not an approximation: no batch's rows ever need reconciling against
// another's. fcs_compute.R (subworkflows/ENRICHMENT/local/src/fcs_compute.R) is
// the extracted, batchable entry point for this; 12.FCS_general_report.Rmd only
// renders the merged result (evidence_score's percentile-rank step and the GMT
// description join both need the FULL merged/cross-database table, so they stay
// in the Rmd, downstream of this).
process FCS_COMPUTE_BATCHED {
    tag "$batchID (${batchSize} GMTs)"
    label 'process_fcs_batched'

    publishDir path: "${params.outdir}/fcs/batches",
               mode: 'copy', overwrite: true,
               enabled: params.publish_intermediates

    input:
    tuple val(batchID), val(batchSize), path(gmtFiles, stageAs: 'gmts/*')
    path stats_file
    path universe_file
    path perms_file

    output:
    path "fcs_enrich_partial.tsv", emit: partial

    script:
    def local_dir = "${baseDir}/subworkflows/ENRICHMENT/local"
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def fdr_permsum     = params.fcs_fdr_permsum     ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def rscript_cmd = (params.use_singularity || params.use_apptainer) ?
        "/usr/local/bin/_entrypoint.sh Rscript" : "Rscript"
    """
    cp ${local_dir}/src/fcs_enrich.R ${local_dir}/src/fcs_compute.R ${local_dir}/src/percentile_flags.R .
    ${rscript_cmd} fcs_compute.R \
        --stats-file ${stats_file} \
        --universe-file ${universe_file} \
        --gmt-dir gmts \
        --perms-file ${perms_file} \
        --num-g ${num_g} \
        --max-g ${max_g} \
        --fdr-thr ${fdr_thr} \
        --fdr-wilcoxon ${fdr_wilcoxon} \
        --fdr-lachenbruch ${fdr_lachenbruch} \
        --fdr-permsum ${fdr_permsum} \
        --pperm-thr ${pperm_thr} \
        --seed ${params.seed ?: 1998} \
        --output fcs_enrich_partial.tsv
    """
}

process FCS_CONCAT {
    tag "Concatenating FCS batch outputs"

    input:
    path(partial_files, stageAs: "partial_*")

    output:
    path "fcs_enrich_merged.tsv", emit: enrich

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mapfile -t files < <(find . -maxdepth 1 -name "partial_*" ! -name ".*" | sort)
    cat "\${files[0]}" > fcs_enrich_merged.tsv
    for ((i=1; i<\${#files[@]}; i++)); do
        tail -n +2 "\${files[\$i]}" >> fcs_enrich_merged.tsv
    done
    """
}

// GMT-batched replacement for calling fcs_run_all() inline inside the Rmd.
// stats_file/universe_file/perms_file are the same inputs SCORING_FCS_REPORT
// / RER_FCS_REPORT already take; params.gmt_dir is the same directory both
// processes already resolve their own `gmt_dir` render param from.
workflow FCS_COMPUTE {
    take:
    stats_file
    universe_file
    perms_file

    main:
    def batchSize = (params.fcs_batch_size ?: 4) as int
    def counter = 0
    def batches = Channel.fromPath("${params.gmt_dir}/*.gmt")
        .collate(batchSize)
        .map { batch ->
            def idx = ++counter
            def batchID = sprintf('fcs_batch_%03d', idx)
            tuple(batchID, batch.size(), batch)
        }

    // stats_file/universe_file/perms_file are take: params -- each carries
    // exactly one item, but crossing this subworkflow's take: boundary loses
    // any value-channel inference Nextflow might have applied upstream (see
    // caas_permulation.nf's CAAS_PERMS_DISAMBIGUATE_BATCHED fix for the full
    // mechanism). Paired positionally against the many-item batches channel,
    // any one of them would silently truncate FCS_COMPUTE_BATCHED to its
    // first batch once exhausted. .first() makes each reusable/broadcastable;
    // no-op if it was already a value channel.
    def stats_file_bc   = stats_file.first()
    def universe_file_bc = universe_file.first()
    def perms_file_bc    = perms_file.first()

    FCS_COMPUTE_BATCHED(batches, stats_file_bc, universe_file_bc, perms_file_bc)
    FCS_CONCAT(FCS_COMPUTE_BATCHED.out.partial.collect())

    emit:
    enrich_file = FCS_CONCAT.out.enrich
}
