#!/usr/bin/env nextflow
// fcs.nf — FCS (Functional Class Scoring) report processes and the merge of the batched tables.
// PhyloPhere | subworkflows/ENRICHMENT/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  SCORING_FCS_REPORT, RER_FCS_REPORT, FCS_CONCAT:
 *  rank-based, threshold-free gene-set enrichment (Wilcoxon-AUC through
 *  RERconverge::fastwilcoxGMTall, with permulation p and BH correction)
 *  over the GMT files of subworkflows/ENRICHMENT/dat.
 *
 *  The report processes render 11.FCS_general_report.Rmd against a generic stats TSV
 *  (gene, score_<ranking> and flag_<name> columns) and a universe file (the
 *  cleaned background for CAAS, genes without signal floored to 0). RER ranks only
 *  genes with a non-missing statistic, retaining measured zeros.
 *
 *    SCORING_FCS_REPORT : CAAS report (global/top/bottom rankings, cross-module flags)
 *    RER_FCS_REPORT     : RERconverge report (also takes the CAAS stats as annot_file)
 *
 *  FADE has no FCS ranking of its own: its statistic is a maximum over many sites,
 *  which does not support a reliable standalone significance test. It enters as a
 *  cross-module corroboration flag on the leading edge of the CAAS and RER rankings.
 *
 *  Consumes:  stats TSV, universe file, permutations file, GMT directory
 *  Produces:  11.FCS_*.html, fcs_results/ (fcs_all_results.tsv, fcs_leading_edge.tsv,
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
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**',
               saveAs: { name -> name.replaceFirst(/^fcs_results\//, '') }

    input:
    path fcs_stats
    path universe
    path enrich_file
    path report_metadata     // fcs_report_metadata.json from FCS_CONCAT; staged beside enrich_file, which the report looks it up by

    output:
    path "11.FCS_scoring_${params.traitname ?: 'unknown_trait'}.html", emit: report
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
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    def enrich_file_arg = (enrich_file.name =~ /^NO_/) ? 'NULL' : "'${enrich_file}'"
    def report_sources = ['11.FCS_general_report.Rmd','src/term_blocks.R','src/fcs_selection.R','src/evidence_flags.R','src/report_presentation.R'].collect { file("${local_dir}/${it}").text }.join('\n')
    def report_hash = java.security.MessageDigest.getInstance('SHA-256').digest(report_sources.getBytes('UTF-8')).encodeHex().toString()
    def render = """
        # FCS reporting source SHA256: ${report_hash}
        rmarkdown::render(
            '11.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = 'Scoring_FCS_${traitname}',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                pperm_thr     = ${pperm_thr},
                analytic_levels = '${params.fcs_analytic_levels}',
                gene_stat = '${params.gene_stat ?: "max"}',
                cluster_threshold  = ${params.fcs_cluster_threshold},
                top_n         = ${top_n},
                traitname     = '${traitname}',
                enrich_file   = ${enrich_file_arg}
            ),
            output_file = '11.FCS_scoring_${traitname}.html'
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
               mode: 'copy', overwrite: true, pattern: 'fcs_results/**',
               saveAs: { name -> name.replaceFirst(/^fcs_results\//, '') }

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
    def pperm_thr = params.fcs_pperm_thr
    def top_n     = params.fcs_top_n
    def enrich_file_arg = (enrich_file.name =~ /^NO_/) ? 'NULL' : "'${enrich_file}'"
    def report_sources = ['11.FCS_general_report.Rmd','src/term_blocks.R','src/fcs_selection.R','src/evidence_flags.R','src/report_presentation.R'].collect { file("${local_dir}/${it}").text }.join('\n')
    def report_hash = java.security.MessageDigest.getInstance('SHA-256').digest(report_sources.getBytes('UTF-8')).encodeHex().toString()
    def render = """
        # FCS reporting source SHA256: ${report_hash}
        rmarkdown::render(
            '11.FCS_general_report.Rmd',
            params = list(
                stats_file    = '${fcs_stats}',
                universe_file = '${universe}',
                gmt_dir       = '${gmt_dir}',
                project_name  = '${report_label}',
                score_source  = 'rer',
                num_g         = ${num_g},
                max_g         = ${max_g},
                fdr_thr       = ${fdr_thr},
                pperm_thr     = ${pperm_thr},
                analytic_levels = '${params.fcs_analytic_levels}',
                gene_stat = '${params.gene_stat ?: "max"}',
                cluster_threshold  = ${params.fcs_cluster_threshold},
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

process FCS_CONCAT {
    tag "Concatenating FCS batch outputs"

    input:
    path(partial_files, stageAs: "partial_*")
    path(calibration_metadata, stageAs: 'calibration_metadata.json')   // [] when not calibrated from SCORING

    output:
    path "fcs_enrich_merged.tsv", emit: enrich
    path "fcs_report_metadata.json", optional: true, emit: report_metadata

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

    # The report restores these settings from this record: the parameters that FCS_OBSERVED and FCS_CLOSE used,
    # which set sig_perm, q.perm and redlich_level of the table. Calibrated runs add the gene_stat that the null
    # check verified upstream (calibration/metadata.json); the report confirms gene_stat=max with it.
    gs=""
    if [ -s calibration_metadata.json ]; then
        gs=\$(sed -nE 's/.*"gene_stat": *"([^"]+)".*/\\1/p' calibration_metadata.json)
        [ -n "\$gs" ] || { echo "FCS_CONCAT: no gene_stat in calibration metadata" >&2; exit 1; }
    fi
    {
        printf '{\\n'
        [ -z "\$gs" ] || printf '  "gene_stat": "%s",\\n' "\$gs"
        printf '  "num_g": %s,\\n  "max_g": %s,\\n  "fdr_thr": %s,\\n  "pperm_thr": %s,\\n  "analytic_levels": "%s"\\n}\\n' \\
            ${params.fcs_min_genes} ${params.fcs_max_genes ?: 0} ${params.fcs_fdr} ${params.fcs_pperm_thr} "${params.fcs_analytic_levels}"
    } > fcs_report_metadata.json
    """
}
