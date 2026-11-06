#!/usr/bin/env nextflow
/*
 * workflow_map.nf
 *
 * Functions that build the final PhyloPhere workflow-map HTML (workflow_map.html).
 * Included by main.nf; the workflow.onComplete hook calls buildCtx and buildWorkflowMapHtml.
 * Written as plain functions instead of a class in lib/ so that the strict linter resolves them.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

    // ── HTML utilities ─────────────────────────────────────────────────────────

    String htmlEscape(value) {
        def s = value == null ? '' : value.toString()
        s = s.replace('&', '&amp;')
             .replace('<', '&lt;')
             .replace('>', '&gt;')
             .replace('"', '&quot;')
             .replace("'", '&#39;')
        return s
    }

    String normalizePath(String p) {
        if (p == null) return null
        return p.toString().replace('\\', '/')
    }

    String imageDataUriIfExists(String path) {
        if (!path) return null
        def f = new File(path)
        if (!f.exists()) return null

        def lower = f.name.toLowerCase()
        def mime = lower.endsWith('.png') ? 'image/png'
                 : (lower.endsWith('.jpg') || lower.endsWith('.jpeg')) ? 'image/jpeg'
                 : 'application/octet-stream'

        def encoded = java.util.Base64.getEncoder().encodeToString(f.bytes)
        return "data:${mime};base64,${encoded}"
    }

    String relativeFromOutdir(String outdir, String targetPath) {
        if (!outdir || !targetPath) return targetPath
        def outNorm = normalizePath(new File(outdir).absolutePath)
        def tgtNorm = normalizePath(new File(targetPath).absolutePath)

        if (tgtNorm == outNorm) return '.'

        def outWithSlash = outNorm.endsWith('/') ? outNorm : "${outNorm}/"
        if (tgtNorm.startsWith(outWithSlash)) {
            return tgtNorm.substring(outWithSlash.length())
        }

        // Fallback: if target is outside outdir, keep absolute path.
        return targetPath
    }

    // ── Evidence of a run ──────────────────────────────────────────────────────

    // Html report files may carry a trait suffix (e.g. 14.Comparison_report_<trait>.html).
    // An exact match wins; otherwise every file of the same directory with the same prefix.
    List<String> htmlMatches(File file) {
        if (file.isFile()) return [file.absolutePath]
        def parent = file.parentFile
        if (!parent?.isDirectory()) return []
        def base = file.name.replace('.html', '')
        return (parent.listFiles() ?: []).findAll { it.name.startsWith(base) && it.name.endsWith('.html') }
            .collect { it.absolutePath }.sort()
    }

    // True when the target is a file, or a directory holding a file, modified at or after sinceMs.
    // Publishing with mode 'copy' rewrites the file on every run, cached or not, so its mtime
    // tells whether this run produced it (directories alone can be left over from an earlier run).
    boolean hasFileSince(File target, long sinceMs) {
        if (!target.exists()) return false
        if (target.isFile()) return target.lastModified() >= sinceMs
        return (target.listFiles() ?: []).any { hasFileSince(it, sinceMs) }
    }

    boolean stageRan(String outdir, List<String> evidence, long sinceMs) {
        return evidence.any { rel ->
            def target = new File(outdir, rel)
            rel.endsWith('.html')
                ? htmlMatches(target).any { new File(it).lastModified() >= sinceMs }
                : hasFileSince(target, sinceMs)
        }
    }

    // ── Stage card ─────────────────────────────────────────────────────────────

    String stageCard(String outdir, Map st) {
        def ran         = st.ran as boolean
        def nodeColor   = ran ? st.color : '#B8B8B8'
        def borderColor = ran ? '#1F2937' : '#9CA3AF'

        def filesDirs = (st.filesDirs ?: []) as List

        def htmlCandidates = (st.htmlCandidates ?: []) as List
        def htmlRows = htmlCandidates ? htmlCandidates.collect { ht ->
            def matches = htmlMatches(new File(ht.toString()))
            if (!matches) {
                def href = relativeFromOutdir(outdir, ht.toString())
                return "<div class=\"link-row\">→ html: <a href=\"${href}\">${htmlEscape(href)}</a> <span class=\"status missing\">MISSING</span></div>"
            }
            return matches.collect { mp ->
                def href = relativeFromOutdir(outdir, mp)
                "<div class=\"link-row\">→ html: <a href=\"${href}\">${htmlEscape(href)}</a> <span class=\"status ok\">available</span></div>"
            }.join('\n')
        }.join('\n') : ''

        def filesRows = filesDirs ? filesDirs.collect { fd ->
            def exists = new File(fd.toString()).exists()
            def label  = exists ? 'available' : 'MISSING'
            def href   = relativeFromOutdir(outdir, fd.toString())
            def text   = relativeFromOutdir(outdir, fd.toString())
            return "<div class=\"link-row\">→ files: <a href=\"${href}\">${htmlEscape(text)}</a> <span class=\"status ${exists ? 'ok' : 'missing'}\">${label}</span></div>"
        }.join('\n') : '<div class="link-row">→ files: <span class="status missing">MISSING</span></div>'

        return """
    <div class=\"stage\" style=\"border-color:${borderColor};\">
      <div class=\"stage-head\" style=\"background:${nodeColor};\">${htmlEscape(st.name)}</div>
      <div class=\"stage-body\">
        <div><span class=\"pill ${ran ? 'pill-ran' : 'pill-skip'}\">${ran ? 'ran' : 'not run'}</span></div>
        ${filesRows}
        ${htmlRows}
      </div>
    </div>
    """.stripIndent()
    }

    // ── Full HTML page ─────────────────────────────────────────────────────────

    String buildWorkflowMapHtml(Map ctx) {
        def outdir     = ctx.outdir
        def projectDir = ctx.projectDir ?: outdir
        def logoDataUri = imageDataUriIfExists("${projectDir}/res/logo.png")

        def colors = [
            reporting : '#7C3AED',   // reports / Enrichment
            prepost   : '#0EA5E9',   // pre/post-processing / characterization
            processes : '#F97316'    // CT / RER / disambiguation / accumulation / selection tools
        ]

        // One row per stage: the card links (filesDirs, htmlCandidates). Run status comes from
        // stageEvidence(), the single table of evidence paths.
        def stages = [
            [ id: 'prune',       name: 'Data pruning',                    type: 'prepost',
              filesDirs: ["${outdir}/data_exploration/0.Data-pruning"],
              htmlCandidates: ["${outdir}/html_reports/2.Phenotype_exploration_pruned.html"] ],

            [ id: 'dataset_rep', name: 'Dataset reporting',               type: 'reporting',
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/html_reports/1.Dataset_exploration.html"] ],

            [ id: 'pheno_rep',   name: 'Phenotype reporting',             type: 'reporting',
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/html_reports/2.Phenotype_exploration_complete.html"] ],

            // directories.R (TRAIT_ANALYSIS) creates 1.Traitfiles and 2.Permulation_traitfiles under 2.CT.
            [ id: 'contrast',    name: 'Contrast selection',              type: 'prepost',
              filesDirs: ["${outdir}/data_exploration/2.CT",
                          "${outdir}/data_exploration/2.CT/1.Traitfiles",
                          "${outdir}/data_exploration/2.CT/2.Permulation_traitfiles"],
              htmlCandidates: ["${outdir}/html_reports/3.CI-composition.html",
                                "${outdir}/html_reports/4.Independent_contrasts.html"] ],

            // discovery/ and resample/ are published only with --publish_intermediates
            // (conf/ct.config); caastools/ is always written and drives the "ran" status.
            [ id: 'ct',          name: 'CT (convergence)',                type: 'processes',
              filesDirs: ["${outdir}/caastools", "${outdir}/discovery",
                          "${outdir}/resample"],
              htmlCandidates: [] ],

            // CAAS_core publishes the permutation null (caas_permulation/) and, with
            // --strict_conservation_output, the strict analysis under strict/ (prefix /strict).
            [ id: 'caas_null',   name: 'CAAS null (permulation)',         type: 'processes',
              filesDirs: ["${outdir}/caas_permulation"],
              htmlCandidates: [] ],

            [ id: 'strict',      name: 'Strict conservation output',      type: 'processes',
              filesDirs: ["${outdir}/strict"],
              htmlCandidates: [] ],

            // 7.CAAS_pattern_annotation.Rmd (CAAS_META_CAAS_REPORT) publishes under meta_caas/.
            [ id: 'ct_signif',   name: 'CAAS Pattern Annotation',         type: 'reporting',
              filesDirs: ["${outdir}/meta_caas",
                          "${outdir}/meta_caas/meta_caas"],
              htmlCandidates: ["${outdir}/html_reports/7.CAAS_pattern_annotation.html"] ],

            [ id: 'ct_disambig', name: 'CT disambiguation (convergence)', type: 'processes',
              filesDirs: ["${outdir}/ct_disambiguation"],
              htmlCandidates: [] ],

            [ id: 'asr_diagnostics', name: 'ASR Diagnostics',              type: 'reporting',
              filesDirs: ["${outdir}/asr_diagnostics"],
              htmlCandidates: ["${outdir}/html_reports/9.ASR_diagnostics.html"] ],

            [ id: 'ct_postproc', name: 'CT post-processing',              type: 'prepost',
              filesDirs: ["${outdir}/postproc",
                          "${outdir}/postproc/preprocessed",
                          "${outdir}/postproc/gene_filtering",
                          "${outdir}/postproc/cleaned_backgrounds"],
              htmlCandidates: ["${outdir}/html_reports/8.Characterization_report.html",
                                "${outdir}/html_reports/8.1.Removed_genes_report.html"] ],

            [ id: 'ct_acc',      name: 'POSENRICH level 1: accumulation',   type: 'processes',
              filesDirs: ["${outdir}/posenrich/accumulation",
                          "${outdir}/posenrich/accumulation/top",
                          "${outdir}/posenrich/accumulation/bottom",
                          "${outdir}/posenrich/accumulation/all"],
              htmlCandidates: [] ],

            [ id: 'vep',         name: 'VEP characterization',            type: 'prepost',
              filesDirs: ["${outdir}/vep"],
              htmlCandidates: [] ],

            // rerconverge/ also holds the per-module AMI, so the RER evidence is rer_results only.
            [ id: 'rer',         name: 'RERconverge (RER)',               type: 'processes',
              filesDirs: ["${outdir}/rerconverge",
                          "${outdir}/rerconverge/rer_results",
                          "${outdir}/rerconverge/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/5.RERconverge_report.html"] ],

            // FADE report files are the evidence: selection/fade/top and bottom also hold AMI output.
            [ id: 'fade',        name: 'FADE (selection)',                type: 'processes',
              filesDirs: ["${outdir}/selection/fade",
                          "${outdir}/selection/fade/top",
                          "${outdir}/selection/fade/bottom",
                          "${outdir}/selection/fade/top/gene_lists",
                          "${outdir}/selection/fade/bottom/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/6.FADE_report_top.html",
                               "${outdir}/html_reports/6.FADE_report_bottom.html"] ],

            // scoring_compute.R writes the ranked slices flat as gene_lists/slice_*.tsv.
            [ id: 'scoring',     name: 'CAAS Scoring',                    type: 'reporting',
              filesDirs: ["${outdir}/scoring",
                          "${outdir}/scoring/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/10.Scoring_report.html"] ],

            // CAAS_SIGNIFICANCE_REPORT (15.CAAS_significance_report.Rmd) runs after scoring:
            // it joins the meta_caas table of ct_signif with position_scores.tsv (p.emp, p.adj_bh).
            [ id: 'ct_signif_sig', name: 'CT significance (post-scoring)', type: 'reporting',
              filesDirs: ["${outdir}/meta_caas/significance"],
              htmlCandidates: ["${outdir}/html_reports/15.CAAS_significance_report.html"] ],

            // ENRICHMENT (workflows/enrichment.nf) runs when --scoring and --enrichment are set.
            // FCS covers the CAAS ranking ('fcs') and the RER ranking ('scoring/rer'); FADE has
            // no FCS ranking. The DOMINO-based AMI run (--scoring_ami) publishes to ami/*, and
            // RER and FADE can also publish their own AMI report under their output trees.
            [ id: 'fcs',         name: 'Functional enrichment (FCS)',     type: 'reporting',
              filesDirs: ["${outdir}/fcs",
                          "${outdir}/scoring/rer"],
              htmlCandidates: ["${outdir}/html_reports/11.FCS_scoring.html",
                               "${outdir}/html_reports/11.FCS_rer.html"] ],

            [ id: 'ami',         name: 'AMI (DOMINO active modules)',     type: 'reporting',
              filesDirs: ["${outdir}/ami",
                          "${outdir}/ami/ami_results",
                          "${outdir}/ami/ami_summary",
                          "${outdir}/ami/ami_plots",
                          "${outdir}/ami/ami_networks",
                          // Per-module AMI, run from the RER and FADE workflow files.
                          "${outdir}/rerconverge/ami",
                          "${outdir}/selection/fade/top/ami",
                          "${outdir}/selection/fade/bottom/ami"],
              htmlCandidates: ["${outdir}/html_reports/12.AMI_analysis.html",
                               "${outdir}/html_reports/12.AMI_rer.html",
                               "${outdir}/html_reports/12.AMI_fade.html"] ],

            [ id: 'posenrich', name: 'POSENRICH registry and level 2',   type: 'reporting',
              filesDirs: ["${outdir}/posenrich/registry",
                          "${outdir}/posenrich/characterization",
                          "${outdir}/posenrich/universe",
                          "${outdir}/posenrich/within_gene",
                          "${outdir}/posenrich/proximity"],
              htmlCandidates: ["${outdir}/html_reports/13.POSENRICH_report.html"] ],

            [ id: 'compare',     name: 'Cross-module comparison',         type: 'reporting',
              filesDirs: ["${outdir}/compare",
                          "${outdir}/compare/compare_results"],
              htmlCandidates: ["${outdir}/html_reports/14.Comparison_report.html"] ]

        ].collect { st -> st + [color: colors[st.type], ran: ctx.ran[st.id] ?: false] }

        def chainIds = ['prune','dataset_rep','pheno_rep','contrast','ct','caas_null','strict','ct_signif',
                        'ct_disambig','asr_diagnostics','ct_postproc','rer','fade',
                        'scoring','ct_signif_sig','vep','ct_acc','posenrich','fcs','ami','compare']
        def rows = []
        chainIds.eachWithIndex { sid, idx ->
            def st = stages.find { it.id == sid }
            rows << stageCard(outdir, st)
            if (idx < chainIds.size() - 1) rows << '<div class="arrow">↓</div>'
        }

        def configSources = [
            "${projectDir}/nextflow.config",
            "${projectDir}/conf/resources.config",
            "${projectDir}/conf/common.config",
            "${projectDir}/conf/ct.config",
            "${projectDir}/conf/ct_disambiguation.config",
            "${projectDir}/conf/ct_postproc.config",
            "${projectDir}/conf/enrichment.config",
            "${projectDir}/conf/posenrich.config",
            "${projectDir}/conf/vep.config",
            "${projectDir}/conf/rerconverge.config",
            "${projectDir}/conf/fade.config",
            "${projectDir}/conf/scoring.config"
        ]

        return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>PhyloPhere workflow map</title>
  <style>
    body { font-family: Inter, Arial, sans-serif; margin: 18px; color: #111827; background: #FAFAFA; }
    h1 { margin: 0 0 6px 0; }
    .subtitle { color:#4B5563; margin-bottom: 14px; }
    .hero { margin-bottom: 16px; }
    .logo-wrap { margin: 0 0 18px 0; }
    .hero-logo { display:block; max-width: 240px; width: 100%; height: auto; object-fit: contain; }
    .meta { font-size: 13px; color:#374151; background:#F3F4F6; border:1px solid #E5E7EB; padding:10px; border-radius:8px; }
    .grid { margin-top: 14px; display:flex; flex-direction:column; gap:8px; }
    .arrow { text-align:center; color:#6B7280; font-size: 18px; }
    .stage { border: 2px solid #D1D5DB; border-radius: 10px; background: #fff; overflow: hidden; }
    .stage-head { color: #fff; font-weight: 700; padding: 8px 10px; }
    .stage-body { padding: 9px 10px 10px; font-size: 13px; }
    .link-row { margin-top: 6px; }
    .pill { display:inline-block; padding:2px 8px; border-radius:999px; font-size:12px; font-weight:600; }
    .pill-ran  { background:#DCFCE7; color:#065F46; }
    .pill-skip { background:#E5E7EB; color:#374151; }
    .status { margin-left:8px; font-weight:700; }
    .ok      { color:#047857; }
    .missing { color:#B91C1C; }
    a { color:#1D4ED8; text-decoration:none; }
    a:hover { text-decoration:underline; }
    .legend { margin-top: 16px; padding: 10px; border:1px solid #E5E7EB; border-radius:8px; background:#fff; }
    .sw { display:inline-block; width:14px; height:14px; border-radius:3px; margin-right:6px; vertical-align:middle; }
    .footer { margin-top:14px; font-size:12px; color:#6B7280; }
    code { background:#F3F4F6; padding:1px 5px; border-radius:4px; }
    ul { margin: 8px 0 0 20px; }
    @media (max-width: 780px) {
      .hero-logo { max-width: 180px; }
    }
  </style>
</head>
<body>
  <div class="hero">
    <div class="logo-wrap">
      ${logoDataUri ? "<img class=\"hero-logo\" src=\"${logoDataUri}\" alt=\"PhyloPhere logo\" />" : ""}
    </div>
    <h1>PhyloPhere workflow map</h1>
    <div class="subtitle">Chain of stages in pipeline order, from pruning and reporting to enrichment (always shown). Gray = not run, colored = run.</div>
  </div>

  <div class="meta">
    <div><b>Run directory:</b> <code>${htmlEscape(ctx.launchDir)}</code></div>
    <div><b>Project directory:</b> <code>${htmlEscape(ctx.projectDir)}</code></div>
    <div><b>Outdir:</b> <code>${htmlEscape(outdir)}</code></div>
    <div><b>Profile:</b> <code>${htmlEscape(ctx.profile)}</code></div>
    <div><b>Run name:</b> <code>${htmlEscape(ctx.runName)}</code></div>
    <div><b>Session ID:</b> <code>${htmlEscape(ctx.sessionId)}</code></div>
    <div><b>Command line:</b> <code>${htmlEscape(ctx.commandLine)}</code></div>
  </div>

  <div class="grid">
    ${rows.join('\n')}
  </div>

  <div class="legend">
    <b>Color key (process type)</b><br/>
    <span class="sw" style="background:${colors.reporting};"></span> Reporting-driven (reports, Enrichment, Scoring)
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <span class="sw" style="background:${colors.prepost};"></span> Pre/postprocessing &amp; characterization (VEP)
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <span class="sw" style="background:${colors.processes};"></span> Analysis processes (CT, RER, FADE, disambiguation, POSENRICH level 1)
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <span class="sw" style="background:#B8B8B8;"></span> Not run
  </div>

  <div class="legend">
    <b>Config source + execution context</b>
    <ul>
      ${configSources.collect { "<li><code>${htmlEscape(it)}</code></li>" }.join('\n')}
    </ul>
    <div style="margin-top:8px;">Loaded profile: <code>${htmlEscape(ctx.profile)}</code></div>
    <div>Executed from: <code>${htmlEscape(ctx.launchDir)}</code></div>
  </div>

  <div class="footer">
    Generated automatically as final run artifact: <code>${htmlEscape(outdir)}/workflow_map.html</code>
    <br/>
    Stage status reflects this run: a stage is "ran" when its evidence files were written after the run started.
  </div>
</body>
</html>
""".stripIndent()
    }

    // ── Run status ─────────────────────────────────────────────────────────────

    // Status of every stage. sinceMs is the start of this run (epoch milliseconds); 0 counts any file,
    // which is what a standalone call without workflow metadata gets.
    Map scanRunStatus(String outdir, long sinceMs) {
        return stageEvidence().collectEntries { id, evidence ->
            [(id): stageRan(outdir, evidence, sinceMs)]
        }
    }

    Map stageEvidence() {
        return [
            prune: ['data_exploration/0.Data-pruning'],
            dataset_rep: ['html_reports/1.Dataset_exploration.html'],
            pheno_rep: ['html_reports/2.Phenotype_exploration_complete.html'],
            contrast: ['data_exploration/2.CT/1.Traitfiles'],
            ct: ['caastools'],
            caas_null: ['caas_permulation'],
            strict: ['strict'],
            ct_signif: ['meta_caas/meta_caas'],
            ct_disambig: ['ct_disambiguation'],
            asr_diagnostics: ['asr_diagnostics'],
            ct_postproc: ['postproc'],
            ct_acc: ['posenrich/accumulation'],
            vep: ['vep'],
            rer: ['rerconverge/rer_results'],
            fade: ['html_reports/6.FADE_report_top.html', 'html_reports/6.FADE_report_bottom.html'],
            scoring: ['scoring/gene_lists'],
            ct_signif_sig: ['meta_caas/significance'],
            fcs: ['fcs'],
            ami: ['ami'],
            posenrich: ['posenrich/registry'],
            compare: ['compare', 'html_reports/14.Comparison_report.html']
        ]
    }

    // ── Convenience: build ctx map from the output directory and workflow metadata ─
    // Can be called with just a directory path for standalone operation,
    // or with params + workflow for backward compatibility.
    // One signature (the hook is its only caller): the output directory, the params and the workflow metadata.
    Map buildCtx(baseDir, paramsMap, workflowMeta) {
        def outdir = baseDir?.toString() ?:
                    (paramsMap?.outdir ? paramsMap.outdir.toString() : "${workflowMeta?.projectDir}/out")

        // workflow.start is an OffsetDateTime in Nextflow; without it every file counts.
        def start = workflowMeta?.start
        long sinceMs = start instanceof java.time.OffsetDateTime ? start.toInstant().toEpochMilli() : 0L

        [
            outdir        : outdir,
            launchDir     : workflowMeta?.launchDir?.toString() ?: outdir,
            projectDir    : workflowMeta?.projectDir?.toString() ?: outdir,
            profile       : workflowMeta?.profile ?: 'standalone',
            runName       : workflowMeta?.runName ?: 'dynamic-scan',
            sessionId     : workflowMeta?.sessionId?.toString() ?: java.util.UUID.randomUUID().toString(),
            commandLine   : workflowMeta?.commandLine ?: 'dynamic generation',
            sinceMs       : sinceMs,
            ran           : scanRunStatus(outdir, sinceMs)
        ]
    }
