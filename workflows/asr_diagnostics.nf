#!/usr/bin/env nextflow
// asr_diagnostics.nf — Workflow wrapper of the ASR diagnostics report.
// PhyloPhere | workflows/

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
# File: asr_diagnostics.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  ASR_DIAGNOSTICS: standalone description of the ASR path score of the observed
 *  scoring. It renders the ASR diagnostics report (ASR_DIAGNOSTICS_REPORT), which
 *  describes the distribution of asr_path_score, of its descriptor derived_agreement
 *  and of the ancestral-state posteriors of the Voronoi domains over the scored rows
 *  of caas_convergence_master.csv.
 *
 *  The workflow runs beside the CT_POSTPROC filtering and does not feed it. The
 *  posterior threshold (params.ct_disambig_posterior_threshold) is applied upstream,
 *  when the ASR posteriors are read; the report displays it and does not vary it.
 *
 *  The ct_disambiguation/ directory comes from the upstream channel or, standalone,
 *  from --disambiguation_dir (a directory or a .tar.gz / .tgz archive).
 *
 *  Consumes:  ct_disambiguation/ directory (caas_convergence_master.csv)
 *  Produces:  report (9.ASR_diagnostics.html), tables (tsv/**), plots (plots/**),
 *             published in asr_diagnostics/ and html_reports/
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { ASR_DIAGNOSTICS_REPORT } from '../subworkflows/ASR_DIAGNOSTICS/asr_diagnostics'

// ── Archive input ────────────────────────────────────────────────────────────

// Unpack a .tar.gz archive into extracted/, dropping its top-level directory.
// Used only when --disambiguation_dir points to an archive.
process EXTRACT_DISAMBIG_DIR {
    input:  path tarball
    output: path "extracted", emit: dir
    script: "mkdir -p extracted && tar -xzf '${tarball}' --strip-components=1 -C extracted/"
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow ASR_DIAGNOSTICS {
    take:
        disambiguation_dir_channel    // ct_disambiguation/ directory of the observed scoring (holds caas_convergence_master.csv)

    main:
        // The disambiguation directory: the upstream channel, else --disambiguation_dir.
        def disambig_dir_ch

        if (disambiguation_dir_channel) {
            log.info "📥 [asr_diagnostics] Using the observed scoring output from upstream"
            disambig_dir_ch = disambiguation_dir_channel
        } else {
            assert params.disambiguation_dir : \
                "[asr_diagnostics] Requires --ct_disambiguation upstream or --disambiguation_dir (path to ct_disambiguation/ directory)"
            def d = file(params.disambiguation_dir)
            assert d.exists() : "[asr_diagnostics] disambiguation_dir not found: ${params.disambiguation_dir}"
            def dStr = params.disambiguation_dir as String
            if (dStr.endsWith('.tar.gz') || dStr.endsWith('.tgz')) {
                disambig_dir_ch = EXTRACT_DISAMBIG_DIR(Channel.value(d)).dir.collect().map { it[0] }
            } else {
                disambig_dir_ch = Channel.value(d)
            }
        }

        def threshold_ch = Channel.value(params.ct_disambig_posterior_threshold)

        log.info "🔬 [asr_diagnostics] Posterior threshold (params.ct_disambig_posterior_threshold): ${params.ct_disambig_posterior_threshold}"
        log.info "📊 [asr_diagnostics] Outputs → ${params.outdir}/asr_diagnostics/"

        diagnostics_output = ASR_DIAGNOSTICS_REPORT(disambig_dir_ch, threshold_ch)

    emit:
        report  = diagnostics_output.report   // 9.ASR_diagnostics.html
        tables  = diagnostics_output.tables   // tsv/**
        plots   = diagnostics_output.plots    // plots/**
}
