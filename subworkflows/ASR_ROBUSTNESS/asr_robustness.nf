#!/usr/bin/env nextflow

/*
#  ASR Robustness: Characterization and reporting (Rmarkdown)
#
#  Consumes the full ct_disambiguation/ output directory produced by
#  CT_DISAMBIGUATION_RUN.  Focal MRCA nodes are auto-derived from the
#  mrca_1_node / mrca_2_node / … columns in caas_convergence_master.csv.
#  The canonical filtering threshold is params.ct_disambig_posterior_threshold.
*/

process ASR_ROBUSTNESS_REPORT {
    tag "asr_robustness"
    label 'process_reporting'
    publishDir path: "${params.outdir}/asr_robustness", mode: 'copy', overwrite: true, pattern: '{tsv/**,plots/**}'
    publishDir path: "${params.outdir}/HTML_reports",   mode: 'copy', overwrite: true, pattern: '*.html'

    input:
    path disambiguation_dir   // full ct_disambiguation/ output directory
    val  posterior_threshold  // params.ct_disambig_posterior_threshold

    output:
    path "*.html",  emit: report
    path "tsv/**",  emit: tables,  optional: true
    path "plots/**", emit: plots,  optional: true

    script:
    def local_dir         = "${baseDir}/subworkflows/ASR_ROBUSTNESS/local"
    def disambig_dir_str  = disambiguation_dir.toString()
    def threshold_str     = posterior_threshold.toString()
    def outdir_str        = "${params.outdir}/asr_robustness"

    if (params.use_singularity || params.use_apptainer) {
        """
        cp ${local_dir}/ASR_robustness.Rmd .

        /usr/local/bin/_entrypoint.sh Rscript -e "
            rmarkdown::render(
                'ASR_robustness.Rmd',
                params = list(
                    disambig_dir        = '${disambig_dir_str}',
                    posterior_threshold = ${threshold_str},
                    output_dir          = '.'
                ),
                output_file = 'ASR_robustness.html'
            )
        "
        """
    } else {
        """
        cp ${local_dir}/ASR_robustness.Rmd .

        Rscript -e "
            rmarkdown::render(
                'ASR_robustness.Rmd',
                params = list(
                    disambig_dir        = '${disambig_dir_str}',
                    posterior_threshold = ${threshold_str},
                    output_dir          = '.'
                ),
                output_file = 'ASR_robustness.html'
            )
        "
        """
    }
}
