#!/usr/bin/env nextflow

/*
 * ACCUMULATION_REPORT
 * ───────────────────
 * Render an HTML report from CT_ACCUMULATION randomization outputs.
 * Receives all per-direction per-scheme CSVs (flat, from RANDOMIZE) plus the
 * aggregation CSVs (from AGGREGATE), reconstructs the expected directory
 * structure in the work directory, then calls 10.Accumulation_report.Rmd.
 *
 * Inputs
 * ──────
 *   rand_csvs  : path — all accumulation_{dir}_{scheme}_aggregated_results.csv files (collected)
 *   agg_csvs   : path — accumulation_global.csv (collected)
 *
 * Outputs
 * ───────
 *   report      : path — 10.Accumulation_report.html
 *   summary_tsv : path — accumulation_summary_{trait}.tsv (optional)
 */

process ACCUMULATION_REPORT {
    tag "accumulation_report|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'ignore'

    publishDir path: "${params.outdir}/accumulation/aggregation",
               mode: 'copy', overwrite: true,
               pattern: '*.html'
    publishDir path: "${params.outdir}/accumulation/aggregation",
               mode: 'copy', overwrite: true,
               pattern: '*.tsv'
    publishDir path: "${params.outdir}/html_reports",
               mode: 'copy', overwrite: true,
               pattern: '*.html'

    input:
    path rand_csvs   // all *_aggregated_results.csv files staged flat
    path agg_csvs    // *_global.csv staged flat

    output:
    path "10.Accumulation_report.html",      emit: report
    path "accumulation_summary_*.tsv",    emit: summary_tsv, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/CT_ACCUMULATION/local"
    def outdir    = "${params.outdir}/accumulation/aggregation"
    def traitname = params.traitname ?: 'unknown_trait'
    def pval_thr  = params.accumulation_report_pval_threshold ?: 0.05

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

        # Reconstruct directory structure expected by 10.Accumulation_report.Rmd
        mkdir -p accum_root/top/randomization \\
                 accum_root/bottom/randomization \\
                 accum_root/all/randomization \\
                 accum_root/aggregation

        for f in accumulation_top_*.csv;    do [ -f "\$f" ] && cp "\$f" accum_root/top/randomization/;    done
        for f in accumulation_bottom_*.csv; do [ -f "\$f" ] && cp "\$f" accum_root/bottom/randomization/; done
        for f in accumulation_all_*.csv;    do [ -f "\$f" ] && cp "\$f" accum_root/all/randomization/;    done
        for f in *_global.csv; do [ -f "\$f" ] && cp "\$f" accum_root/aggregation/ || true; done

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '10.Accumulation_report.Rmd',
                params = list(
                    accum_dir      = 'accum_root',
                    traitname      = '${traitname}',
                    pval_threshold = ${pval_thr},
                    output_dir     = '${outdir}'
                ),
                output_file = '10.Accumulation_report.html'
            )
        "
        """
    } else {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

        mkdir -p accum_root/top/randomization \\
                 accum_root/bottom/randomization \\
                 accum_root/all/randomization \\
                 accum_root/aggregation

        for f in accumulation_top_*.csv;    do [ -f "\$f" ] && cp "\$f" accum_root/top/randomization/;    done
        for f in accumulation_bottom_*.csv; do [ -f "\$f" ] && cp "\$f" accum_root/bottom/randomization/; done
        for f in accumulation_all_*.csv;    do [ -f "\$f" ] && cp "\$f" accum_root/all/randomization/;    done
        for f in *_global.csv; do [ -f "\$f" ] && cp "\$f" accum_root/aggregation/ || true; done

        Rscript -e "
            rmarkdown::render(
                '10.Accumulation_report.Rmd',
                params = list(
                    accum_dir      = 'accum_root',
                    traitname      = '${traitname}',
                    pval_threshold = ${pval_thr},
                    output_dir     = '${outdir}'
                ),
                output_file = '10.Accumulation_report.html'
            )
        "
        """
    }
}
