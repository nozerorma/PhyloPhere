#!/usr/bin/env nextflow

// accum_report.nf — Render the HTML accumulation report.
// PhyloPhere | subworkflows/CT_ACCUMULATION/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ACCUMULATION_REPORT: renders 10.Accumulation_report.Rmd from the accumulation outputs.
 *
 *  The randomization CSVs arrive staged flat; the script copies them into the directory
 *  layout the Rmd reads (accum_root/{top,bottom,all}/randomization and accum_root/aggregation)
 *  before rendering. The two script branches differ only in the container entrypoint.
 *
 *  Consumes:  all accumulation_<direction>_<scheme>_aggregated_results.csv files (collected)
 *             and the *_global.csv files of the aggregation (collected)
 *  Produces:  10.Accumulation_report.html, accumulation_summary_*.tsv (optional)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Report rendering ─────────────────────────────────────────────────────────

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
    def fdr_thr   = params.accumulation_fdr ?: 0.1
    def rand_type = params.accumulation_randomization_type ?: 'cons_decile'

    if (params.use_singularity || params.use_apptainer) {
        """
        cp -R ${local_dir}/* .
        find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

        # Rebuild the directory layout that 10.Accumulation_report.Rmd reads
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
                    fdr_threshold  = ${fdr_thr},
                    randomization_type = '${rand_type}',
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
                    fdr_threshold  = ${fdr_thr},
                    randomization_type = '${rand_type}',
                    output_dir     = '${outdir}'
                ),
                output_file = '10.Accumulation_report.html'
            )
        "
        """
    }
}
