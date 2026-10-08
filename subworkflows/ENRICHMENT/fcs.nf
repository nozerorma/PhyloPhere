#!/usr/bin/env nextflow
// fcs.nf — FCS (Functional Class Scoring) report and batched-statistics processes.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_FCS_REPORT, RER_FCS_REPORT, FCS_COMPUTE_BATCHED, FCS_CONCAT, FCS_COMPUTE:
 *  rank-based, threshold-free gene-set enrichment (Wilcoxon-AUC through
 *  RERconverge::fastwilcoxGMTall, plus the Lachenbruch two-part test of
 *  fcs_enrich.R) over the GMT files of subworkflows/ENRICHMENT/dat.
 *
 *  The report processes render 12.FCS_general_report.Rmd against a generic stats TSV
 *  (gene, score_<ranking> and flag_<name> columns) and a universe file (the
 *  cleaned background, genes without signal floored to 0).
 *
 *    SCORING_FCS_REPORT : CAAS report (global/top/bottom rankings, cross-module flags)
 *    RER_FCS_REPORT     : RERconverge report (also takes the CAAS stats as annot_file)
 *
 *  FADE and Accumulation have no FCS ranking of their own: the FADE statistic is a
 *  maximum over many sites and Accumulation has no permulation null, so neither
 *  supports a reliable standalone significance test. They enter as cross-module
 *  corroboration flags on the leading edge of the CAAS and RER rankings.
 *
 *  Consumes:  stats TSV, universe file, permutations file, GMT directory
 *  Produces:  12.FCS_*.html, fcs_results/ (fcs_all_results.tsv, fcs_leading_edge.tsv,
 *             fcs_leading_edge_composition.tsv), fcs_enrich_merged.tsv
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── SCORING_FCS_REPORT: CAAS FCS report ──────────────────────────────────────
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
    def gmt_dir   = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    // Published gene_lists/ of the scoring module. This is the CAAS report, so the
    // score_top/score_bottom columns are the CAAS rankings (see the gene_lists_dir
    // parameter of 12.FCS_general_report.Rmd); a NO_* sentinel means none was given.
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
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${traitname}',
                perms_file    = '${perms_file}',
                gene_lists_dir = ${gene_lists_arg},
                enrich_file   = ${enrich_file_arg}
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

// ── RER_FCS_REPORT: RERconverge FCS report ───────────────────────────────────
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
    def gmt_dir   = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def num_g     = params.fcs_min_genes
    def max_g     = params.fcs_max_genes ?: 0
    def fdr_thr   = params.fcs_fdr
    def fdr_wilcoxon    = params.fcs_fdr_wilcoxon    ?: params.fcs_fdr
    def fdr_lachenbruch = params.fcs_fdr_lachenbruch ?: params.fcs_fdr
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
                pperm_thr     = ${pperm_thr},
                top_n         = ${top_n},
                traitname     = '${params.traitname ?: "trait"}',
                perms_file    = '${perms_file}',
                annot_file    = '${annot_file}',
                enrich_file   = ${enrich_file_arg}
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

// ── Batched FCS statistics ───────────────────────────────────────────────────
// fcs_run_all() over every GMT database is the expensive step behind both report
// processes (Wilcoxon-AUC and Lachenbruch for every
// score_<ranking> column against every GMT). The BH correction of fcs_enrich.R is
// scoped per database, so the GMT set can be split across independent tasks and the
// partial tables row-concatenated without any reconciliation: the result is exact.
// fcs_compute.R (subworkflows/ENRICHMENT/local/src/) is the batchable entry point.
// 12.FCS_general_report.Rmd only renders the merged table: the evidence_score
// percentile ranking and the GMT description join need the full cross-database
// table, so they stay in the Rmd.
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
        --pperm-thr ${pperm_thr} \
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

    # Batches can write their columns in different orders (a zero-hit batch
    # emits the empty-result schema of fcs_run_all()), so rows are merged by column
    # name onto the header of the first batch that has data rows.
    ref="\${files[0]}"
    for f in "\${files[@]}"; do
        if [ "\$(wc -l < "\$f")" -gt 1 ]; then ref="\$f"; break; fi
    done

    awk -F'\\t' -v OFS='\\t' '
        FNR == 1 {
            delete src
            for (i = 1; i <= NF; i++) src[\$i] = i
            if (NR == 1) {
                n = NF
                for (i = 1; i <= NF; i++) { canon[i] = \$i; want[\$i] = 1 }
                print; next
            }
            for (c in src) if (!(c in want)) {
                printf "FCS_CONCAT: column %s in %s not in reference header\\n", c, FILENAME > "/dev/stderr"
                exit 1
            }
            next
        }
        {
            line = ""
            for (i = 1; i <= n; i++) {
                v = (canon[i] in src) ? \$(src[canon[i]]) : "NA"
                line = (i == 1) ? v : line OFS v
            }
            print line
        }
    ' "\$ref" \$(printf '%s\\n' "\${files[@]}" | grep -vxF "\$ref" || true) > fcs_enrich_merged.tsv
    """
}

// ── FCS_COMPUTE: GMT-batched fcs_run_all() ───────────────────────────────────────
// Takes the same stats, universe and permutations files as the report processes and
// reads the GMT files of params.gmt_dir (the directory the reports use as gmt_dir).
workflow FCS_COMPUTE {
    take:
    stats_file
    universe_file
    perms_file

    main:
    def batchSize = (params.fcs_batch_size ?: 4) as int
    def counter = 0
    def gmt_dir_resolved = params.gmt_dir ?: "${baseDir}/subworkflows/ENRICHMENT/dat"
    def batches = Channel.fromPath("${gmt_dir_resolved}/*.gmt")
        .collate(batchSize)
        .map { batch ->
            def idx = ++counter
            def batchID = String.format('fcs_batch_%03d', idx)
            tuple(batchID, batch.size(), batch)
        }

    // The three inputs carry one item each but are queue channels once they cross a
    // take: boundary. A process pairs its input channels positionally and stops when
    // one runs out, so paired against the many-item batches channel each of them would
    // end FCS_COMPUTE_BATCHED after its first batch (the same mechanism as in
    // CAAS_CORE of subworkflows/CT/caas_permulation.nf). .collect().map { it[0] }
    // turns each into a value channel that every batch can read.
    def stats_file_bc    = stats_file.collect().map { it[0] }
    def universe_file_bc = universe_file.collect().map { it[0] }
    def perms_file_bc    = perms_file.collect().map { it[0] }

    FCS_COMPUTE_BATCHED(batches, stats_file_bc, universe_file_bc, perms_file_bc)
    FCS_CONCAT(FCS_COMPUTE_BATCHED.out.partial.collect())

    emit:
    enrich_file = FCS_CONCAT.out.enrich
}
