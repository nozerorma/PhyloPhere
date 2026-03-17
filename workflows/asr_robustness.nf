#!/usr/bin/env nextflow

/*
##
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: asr_robustness.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ASR Robustness Workflow: Standalone diagnostic module that characterises ASR
 *  posterior uncertainty across all genes and nodes, extracts focal-MRCA confidence
 *  metrics, applies the canonical site-retention filter
 *  (params.ct_disambig_posterior_threshold), and generates publication-quality
 *  diagnostic plots and tables.
 *
 *  Design:
 *    - Runs in parallel with CT_POSTPROC; does NOT alter the CT_POSTPROC data path.
 *    - Focal MRCA nodes are auto-derived from mrca_1_node / mrca_2_node / …
 *      columns in caas_convergence_master.csv (no user-supplied node list needed).
 *    - All filtering uses params.ct_disambig_posterior_threshold exclusively.
 *      No threshold is hardcoded here or in the companion Rmd.
 *    - A sensitivity analysis over tau = {0.90, 0.95, 0.99} is produced as a
 *      diagnostic table; it does not affect pipeline filtering decisions.
 *
 *  Core filtering rule (from the spec):
 *    ASR uncertainty is described globally and per gene, but filtering for
 *    convergence inference is based only on the focal MRCA nodes that define
 *    change polarity. A site is retained only if the minimum MAP posterior
 *    across all focal nodes is >= params.ct_disambig_posterior_threshold.
 *    Additional threshold values are evaluated only in sensitivity analyses.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { ASR_ROBUSTNESS_REPORT } from "${baseDir}/subworkflows/ASR_ROBUSTNESS/asr_robustness"

workflow ASR_ROBUSTNESS {
    take:
        disambiguation_dir_channel    // results_dir from CT_DISAMBIGUATION_RUN (full ct_disambiguation/ dir)

    main:
        // Resolve the disambiguation directory from upstream or standalone param
        def disambig_dir_ch

        if (disambiguation_dir_channel) {
            log.info "📥 [asr_robustness] Using disambiguation output from upstream CT_DISAMBIGUATION"
            disambig_dir_ch = disambiguation_dir_channel
        } else {
            assert params.disambiguation_dir : \
                "[asr_robustness] Requires --ct_disambiguation upstream or --disambiguation_dir (path to ct_disambiguation/ directory)"
            def d = file(params.disambiguation_dir)
            assert d.exists() : "[asr_robustness] disambiguation_dir not found: ${params.disambiguation_dir}"
            disambig_dir_ch = Channel.value(d)
        }

        // Require diagnostics JSONL files to be present (ct_disambig_run_diagnostics must be true)
        // This is the default; warn if it may have been disabled.
        if (params.containsKey('ct_disambig_run_diagnostics') && !params.ct_disambig_run_diagnostics) {
            log.warn "[asr_robustness] ct_disambig_run_diagnostics=false — posterior JSONL files will be absent; " +
                     "global/per-gene MAP distributions will be skipped. Enable ct_disambig_run_diagnostics for full output."
        }

        def threshold_ch = Channel.value(params.ct_disambig_posterior_threshold)

        log.info "🔬 [asr_robustness] Posterior threshold (params.ct_disambig_posterior_threshold): ${params.ct_disambig_posterior_threshold}"
        log.info "📊 [asr_robustness] Outputs → ${params.outdir}/asr_robustness/"

        robustness_output = ASR_ROBUSTNESS_REPORT(disambig_dir_ch, threshold_ch)

    emit:
        report  = robustness_output.report   // ASR_robustness.html
        tables  = robustness_output.tables   // tsv/**
        plots   = robustness_output.plots    // plots/**
}
