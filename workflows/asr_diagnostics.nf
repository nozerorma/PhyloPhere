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
 *  over the rows of caas_convergence_master.csv, with validated observed provenance.
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

// Nextflow deep caching does not hash files inside a directory. Only these two
// scientific inputs are consumed by the report; stream them to bound memory.
def asrInputHash(directory) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    ['caas_convergence_master.csv', 'observed_scoring.json'].each { name ->
        def source = directory.resolve(name)
        assert source.isFile() : "[asr_diagnostics] Required observed input not found: ${source}"
        md.update(name.getBytes('UTF-8'))
        md.update((byte)0)
        source.toFile().withInputStream { stream ->
            byte[] buffer = new byte[1048576]
            int n
            while ((n = stream.read(buffer)) != -1) md.update(buffer, 0, n)
        }
        md.update((byte)0)
    }
    return md.digest().encodeHex().toString()
}

// ── Archive input ────────────────────────────────────────────────────────────

// Select the unique observed contract; archives may be flat or nested.
def asrShellQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def extraction_source_hash = java.security.MessageDigest.getInstance('SHA-256')
    .digest(file("${baseDir}/subworkflows/ASR_DIAGNOSTICS/local/src/extract_diagnostics.py").bytes).encodeHex().toString()
process EXTRACT_DISAMBIG_DIR {
    ext extractor_hash: extraction_source_hash
    cache 'deep'
    input: path tarball
    output: path 'extracted', emit: dir
    script:
    def helper = "${baseDir}/subworkflows/ASR_DIAGNOSTICS/local/src/extract_diagnostics.py"
    """
    echo '${task.ext.extractor_hash}' > extractor.sha256
    python3 ${asrShellQuote(helper)} ${asrShellQuote(tarball)} extracted
    """
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
                disambig_dir_ch = EXTRACT_DISAMBIG_DIR(Channel.value(d)).dir
            } else {
                assert d.isDirectory() : "[asr_diagnostics] Expected a directory or .tar.gz/.tgz archive: ${d}"
                disambig_dir_ch = Channel.value(d)
            }
        }

        def threshold_ch = Channel.value(params.ct_disambig_posterior_threshold)

        log.info "🔬 [asr_diagnostics] Current posterior threshold (the report reads the recorded value): ${params.ct_disambig_posterior_threshold}"
        log.info "📊 [asr_diagnostics] Outputs → ${params.outdir}/asr_diagnostics/"

        diagnostics_output = ASR_DIAGNOSTICS_REPORT(
            disambig_dir_ch.map { d -> tuple(asrInputHash(d), d) }, threshold_ch)

    emit:
        report  = diagnostics_output.report   // 9.ASR_diagnostics.html
        tables  = diagnostics_output.tables   // tsv/**
        plots   = diagnostics_output.plots    // plots/**
}
