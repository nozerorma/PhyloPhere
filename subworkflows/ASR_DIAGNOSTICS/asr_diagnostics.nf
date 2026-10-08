#!/usr/bin/env nextflow
// asr_diagnostics.nf — Render the ASR diagnostics report from the observed scoring output.
// PhyloPhere | subworkflows/ASR_DIAGNOSTICS/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ASR_DIAGNOSTICS_REPORT: renders 9.ASR_diagnostics.Rmd, which describes the
 *  distribution of asr_path_score, of its descriptor derived_agreement and of the
 *  ancestral-state posteriors of the Voronoi domains over the scored rows of
 *  caas_convergence_master.csv.
 *
 *  The report only displays params.ct_disambig_posterior_threshold; the threshold is
 *  applied upstream, when the ASR posteriors are read, and is not varied here.
 *
 *  Consumes:  ct_disambiguation/ directory of the observed scoring (CAAS_CORE_OBSERVED, or
 *             CAAS_OBSERVED when a discovery.tab is reused), posterior threshold
 *  Produces:  9.ASR_diagnostics.html, tsv/ (summary tables), plots/ (PNG)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── ASR diagnostics report ─────────────────────────────────────────────────────

process ASR_DIAGNOSTICS_REPORT {
    tag "asr_diagnostics"
    label 'process_reporting'
    label 'error_retry'
    publishDir path: "${params.outdir}/asr_diagnostics", mode: 'copy', overwrite: true, pattern: '{tsv/**,plots/**}'
    publishDir path: "${params.outdir}/html_reports",   mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path disambiguation_dir   // full ct_disambiguation/ output directory
    val  posterior_threshold  // params.ct_disambig_posterior_threshold

    output:
    path "*.html",  emit: report
    path "tsv/**",  emit: tables,  optional: true
    path "plots/**", emit: plots,  optional: true

    script:
    def local_dir         = "${baseDir}/subworkflows/ASR_DIAGNOSTICS/local"
    def disambig_dir_str  = disambiguation_dir.toString()
    def threshold_str     = posterior_threshold.toString()
    def outdir_str        = "${params.outdir}/asr_diagnostics"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/9.ASR_diagnostics.Rmd .

        REPORT_CORES=${task.cpus} /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                '9.ASR_diagnostics.Rmd',
                params = list(
                    disambig_dir        = '${disambig_dir_str}',
                    posterior_threshold = ${threshold_str},
                    output_dir          = '.'
                ),
                output_file = '9.ASR_diagnostics.html'
            )
        "
        """
    } else {
        """
        cp ${local_dir}/9.ASR_diagnostics.Rmd .

        REPORT_CORES=${task.cpus} Rscript -e "
            rmarkdown::render(
                '9.ASR_diagnostics.Rmd',
                params = list(
                    disambig_dir        = '${disambig_dir_str}',
                    posterior_threshold = ${threshold_str},
                    output_dir          = '.'
                ),
                output_file = '9.ASR_diagnostics.html'
            )
        "
        """
    }
}
