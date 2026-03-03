/*
 * WorkflowMap.groovy
 *
 * Helper library for generating the final PhyloPhere workflow-map HTML artifact.
 * Placed in lib/ so it is automatically compiled and available to all .nf files
 * in the project (main.nf's workflow.onComplete and workflows/workflow_map.nf).
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

import java.util.UUID

class WorkflowMap {

    // ── HTML utilities ─────────────────────────────────────────────────────────

    static String htmlEscape(value) {
        def s = value == null ? '' : value.toString()
        s = s.replace('&', '&amp;')
             .replace('<', '&lt;')
             .replace('>', '&gt;')
             .replace('"', '&quot;')
             .replace("'", '&#39;')
        return s
    }

    static String resolveFirstExisting(List<String> candidates) {
        if (!candidates) return null
        for (def c : candidates) {
            if (new File(c).exists()) return c
        }
        return candidates[0]
    }

    static String normalizePath(String p) {
        if (p == null) return null
        return p.toString().replace('\\', '/')
    }

    static String relativeFromOutdir(String outdir, String targetPath) {
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

    // ── Stage card ─────────────────────────────────────────────────────────────

    static String stageCard(String outdir, Map st) {
        def ran         = st.ran as boolean
        def nodeColor   = ran ? st.color : '#B8B8B8'
        def borderColor = ran ? '#1F2937' : '#9CA3AF'

        def filesDirs = (st.filesDirs ?: []) as List
        if (!filesDirs && st.filesDir) filesDirs = [st.filesDir]

        def htmlTarget = resolveFirstExisting(st.htmlCandidates ?: [])
        def htmlExists = htmlTarget ? new File(htmlTarget).exists() : false
        def htmlLabel  = htmlExists ? 'available' : 'MISSING'
        def htmlHref   = htmlTarget ? relativeFromOutdir(outdir, htmlTarget) : '#'
        def htmlText   = htmlTarget ? relativeFromOutdir(outdir, htmlTarget) : '(none)'

        def filesRows = filesDirs ? filesDirs.collect { fd ->
            def exists = new File(fd.toString()).exists()
            def label  = exists ? 'available' : 'MISSING'
            def href   = relativeFromOutdir(outdir, fd.toString())
            def text   = relativeFromOutdir(outdir, fd.toString())
            return "<div class=\"link-row\">\u2192 files: <a href=\"${href}\">${htmlEscape(text)}</a> <span class=\"status ${exists ? 'ok' : 'missing'}\">${label}</span></div>"
        }.join('\n') : '<div class="link-row">\u2192 files: <span class="status missing">MISSING</span></div>'

        return """
    <div class=\"stage\" style=\"border-color:${borderColor};\">
      <div class=\"stage-head\" style=\"background:${nodeColor};\">${htmlEscape(st.name)}</div>
      <div class=\"stage-body\">
        <div><span class=\"pill ${ran ? 'pill-ran' : 'pill-skip'}\">${ran ? 'ran' : 'not run'}</span></div>
        ${filesRows}
        <div class=\"link-row\">\u2192 html: <a href=\"${htmlHref}\">${htmlEscape(htmlText)}</a> <span class=\"status ${htmlExists ? 'ok' : 'missing'}\">${htmlLabel}</span></div>
      </div>
    </div>
    """.stripIndent()
    }

    // ── Full HTML page ─────────────────────────────────────────────────────────

    static String buildWorkflowMapHtml(Map ctx) {
        def outdir     = ctx.outdir
        def projectDir = ctx.projectDir ?: outdir

        def colors = [
            reporting : '#7C3AED',   // reports / ORA / STRING
            prepost   : '#0EA5E9',   // pre/post-processing
            processes : '#F97316'    // CT / RER / disambiguation / accumulation / selection tools
        ]

        def stages = [
            [ id: 'prune',       name: 'Data pruning',                    type: 'prepost',   ran: ctx.prune,
              filesDirs: ["${outdir}/data_exploration/0.Data-pruning"],
              htmlCandidates: ["${outdir}/HTML_reports/0.Data_pruning.html"] ],

            [ id: 'dataset_rep', name: 'Dataset reporting',               type: 'reporting', ran: ctx.datasetReport,
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/HTML_reports/1.Dataset_exploration.html"] ],

            [ id: 'pheno_rep',   name: 'Phenotype reporting',             type: 'reporting', ran: ctx.phenotypeRep,
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/HTML_reports/2.Phenotype_exploration.html"] ],

            [ id: 'contrast',    name: 'Contrast selection',              type: 'prepost',   ran: ctx.contrastSel,
              filesDirs: ["${outdir}/data_exploration/2.CT",
                          "${outdir}/data_exploration/2.CT/1.Traitfiles",
                          "${outdir}/data_exploration/2.CT/2.Bootstrap_traitfiles",
                          "${outdir}/data_exploration/2.CT/3.Tree"],
              htmlCandidates: ["${outdir}/HTML_reports/4.Independent_contrasts.html",
                               "${outdir}/HTML_reports/3.CI-composition.html"] ],

            [ id: 'ct',          name: 'CT (convergence)',                type: 'processes', ran: ctx.ct,
              filesDirs: ["${outdir}/caastools", "${outdir}/discovery",
                          "${outdir}/resample",  "${outdir}/bootstrap"],
              htmlCandidates: ["${outdir}/HTML_reports/ct_overview.html"] ],

            [ id: 'ct_signif',   name: 'CT signification (convergence)',  type: 'reporting', ran: ctx.ctSignif,
              filesDirs: ["${outdir}/signification",
                          "${outdir}/signification/gene_lists",
                          "${outdir}/signification/meta_caas"],
              htmlCandidates: ["${outdir}/HTML_reports/CT_signification.html"] ],

            [ id: 'ct_disambig', name: 'CT disambiguation (convergence)', type: 'processes', ran: ctx.ctDisambig,
              filesDirs: ["${outdir}/ct_disambiguation"],
              htmlCandidates: ["${outdir}/HTML_reports/CT_disambiguation.html"] ],

            [ id: 'ct_postproc', name: 'CT post-processing',              type: 'prepost',   ran: ctx.ctPostproc,
              filesDirs: ["${outdir}/postproc",
                          "${outdir}/postproc/preprocessed",
                          "${outdir}/postproc/gene_filtering",
                          "${outdir}/postproc/cleaned_backgrounds"],
              htmlCandidates: ["${outdir}/HTML_reports/CT_postproc.html"] ],

            [ id: 'ora',         name: 'ORA',                             type: 'reporting', ran: ctx.ora,
              filesDirs: ["${outdir}/ora", "${outdir}/accumulation/ora"],
              htmlCandidates: ["${outdir}/HTML_reports/ORA_general.html",
                               "${outdir}/HTML_reports/ORA_accumulation.html"] ],

            [ id: 'string',      name: 'STRING',                          type: 'reporting', ran: ctx.string,
              filesDirs: ["${outdir}/string", "${outdir}/accumulation/string"],
              htmlCandidates: ["${outdir}/HTML_reports/STRING_general.html",
                               "${outdir}/HTML_reports/STRING_accumulation.html"] ],

            [ id: 'ct_acc',      name: 'CT accumulation (convergence)',   type: 'processes', ran: ctx.ctAccum,
              filesDirs: ["${outdir}/accumulation",
                          "${outdir}/accumulation/aggregation",
                          "${outdir}/accumulation/randomization"],
              htmlCandidates: ["${outdir}/HTML_reports/CT_accumulation.html"] ],

            [ id: 'rer',         name: 'RERconverge (RER)',               type: 'processes', ran: ctx.rer,
              filesDirs: ["${outdir}/rerconverge", "${outdir}/rerconverge/gene_sets"],
              htmlCandidates: ["${outdir}/HTML_reports/RERconverge.html"] ],

            [ id: 'fade',        name: 'FADE (selection)',                type: 'processes', ran: ctx.fade,
              filesDirs: ["${outdir}/selection/fade",
                          "${outdir}/selection/fade/top",
                          "${outdir}/selection/fade/bottom"],
              htmlCandidates: ["${outdir}/selection/fade/top/FADE_top.html",
                               "${outdir}/selection/fade/bottom/FADE_bottom.html"] ],

            [ id: 'molerate',    name: 'Molerate (RER)',                  type: 'processes', ran: ctx.molerate,
              filesDirs: ["${outdir}/selection/molerate",
                          "${outdir}/selection/molerate/top",
                          "${outdir}/selection/molerate/bottom"],
              htmlCandidates: ["${outdir}/selection/molerate/top/MOLERATE_top.html",
                               "${outdir}/selection/molerate/bottom/MOLERATE_bottom.html"] ]

        ].collect { st -> st + [color: colors[st.type]] }

        def chainIds = ['prune','dataset_rep','pheno_rep','contrast','ct','ct_signif',
                        'ct_disambig','ct_postproc','ora','string','ct_acc','rer','fade','molerate']
        def rows = []
        chainIds.eachWithIndex { sid, idx ->
            def st = stages.find { it.id == sid }
            rows << stageCard(outdir, st)
            if (idx < chainIds.size() - 1) rows << '<div class="arrow">\u2193</div>'
        }

        def configSources = [
            "${projectDir}/nextflow.config",
            "${projectDir}/conf/resources.config",
            "${projectDir}/conf/common.config",
            "${projectDir}/conf/ct.config",
            "${projectDir}/conf/rerconverge.config",
            "${projectDir}/conf/ora.config",
            "${projectDir}/conf/ct_postproc.config",
            "${projectDir}/conf/ct_disambiguation.config",
            "${projectDir}/conf/ct_accumulation.config",
            "${projectDir}/conf/fade.config",
            "${projectDir}/conf/molerate.config"
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
    .meta { font-size: 13px; color:#374151; background:#F3F4F6; border:1px solid #E5E7EB; padding:10px; border-radius:8px; }
    .grid { margin-top: 14px; display:flex; flex-direction:column; gap:8px; max-width: 980px; }
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
    .legend { margin-top: 16px; padding: 10px; border:1px solid #E5E7EB; border-radius:8px; background:#fff; max-width:980px; }
    .sw { display:inline-block; width:14px; height:14px; border-radius:3px; margin-right:6px; vertical-align:middle; }
    .footer { margin-top:14px; font-size:12px; color:#6B7280; }
    code { background:#F3F4F6; padding:1px 5px; border-radius:4px; }
    ul { margin: 8px 0 0 20px; }
  </style>
</head>
<body>
  <h1>PhyloPhere workflow map</h1>
  <div class="subtitle">Complete chain from prune/reporting to selection (always shown). Gray = not run, colored = run.</div>

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
    <span class="sw" style="background:${colors.reporting};"></span> Reporting-driven (reports, ORA, STRING)
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <span class="sw" style="background:${colors.prepost};"></span> Pre/postprocessing
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <span class="sw" style="background:${colors.processes};"></span> Processes as such (CT, RER, selection tools, disambiguation, accumulation)
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
    <button onclick="refreshStatus()" style="margin-top:8px; padding:4px 12px; border:1px solid #D1D5DB; border-radius:6px; background:#F9FAFB; cursor:pointer;">Refresh Status</button>
    <span id="refresh-status" style="margin-left:8px; font-style:italic; color:#6B7280;"></span>
  </div>

  <script>
    function refreshStatus() {
      document.getElementById('refresh-status').textContent = 'Refreshing...';
      
      // Get the current HTML file's directory (baseDir)
      const currentPath = window.location.pathname;
      const baseDir = currentPath.substring(0, currentPath.lastIndexOf('/'));
      
      // Define canonical directories for each stage
      const canonicalDirs = {
        'prune': ['data_exploration/0.Data-pruning'],
        'dataset_rep': ['data_exploration'],
        'pheno_rep': ['data_exploration'],
        'contrast': ['data_exploration/2.CT', 'data_exploration/2.CT/1.Traitfiles'],
        'ct': ['caastools', 'discovery'],
        'ct_signif': ['signification', 'signification/gene_lists'],
        'ct_disambig': ['ct_disambiguation'],
        'ct_postproc': ['postproc', 'postproc/preprocessed'],
        'ora': ['ora'],
        'string': ['string'],
        'ct_acc': ['accumulation', 'accumulation/aggregation'],
        'rer': ['rerconverge'],
        'fade': ['selection/fade'],
        'molerate': ['selection/molerate']
      };
      
      // Map stage IDs to display stages
      const stageMapping = {
        'prune': 'prune',
        'dataset_rep': 'dataset_rep',
        'pheno_rep': 'pheno_rep', 
        'contrast': 'contrast',
        'ct': 'ct',
        'ct_signif': 'ct_signif',
        'ct_disambig': 'ct_disambig',
        'ct_postproc': 'ct_postproc',
        'ora': 'ora',
        'string': 'string',
        'ct_acc': 'ct_acc',
        'rer': 'rer',
        'fade': 'fade',
        'molerate': 'molerate'
      };
      
      let checkedStages = 0;
      const totalStages = Object.keys(canonicalDirs).length;
      
      // Check each stage's directories
      Object.entries(canonicalDirs).forEach(([stageId, dirPaths]) => {
        let stageCompleted = false;
        let checkedPaths = 0;
        
        dirPaths.forEach(dirPath => {
          const fullPath = baseDir + '/' + dirPath;
          
          // Try to detect directory existence (limited by browser security)
          fetch(fullPath + '/')
            .then(response => {
              if (response.ok || response.status === 403) {
                stageCompleted = true;
              }
            })
            .catch(() => {
              // Directory doesn't exist or can't access
            })
            .finally(() => {
              checkedPaths++;
              if (checkedPaths === dirPaths.length) {
                updateStageDisplay(stageId, stageCompleted);
                checkedStages++;
                if (checkedStages === totalStages) {
                  document.getElementById('refresh-status').textContent = 'Refreshed at ' + new Date().toLocaleTimeString();
                }
              }
            });
        });
      });
    }
    
    function updateStageDisplay(stageId, completed) {
      const stages = document.querySelectorAll('.stage');
      stages.forEach(stage => {
        const stageHead = stage.querySelector('.stage-head');
        if (stageHead && stageHead.textContent.includes(getStageNameFromId(stageId))) {
          const pill = stage.querySelector('.pill');
          const stageElement = stage;
          
          if (completed) {
            pill.className = 'pill pill-ran';
            pill.textContent = 'ran';
            stageElement.style.borderColor = '#1F2937';
            stageHead.style.background = getStageColor(stageId);
          } else {
            pill.className = 'pill pill-skip';
            pill.textContent = 'not run';
            stageElement.style.borderColor = '#9CA3AF';
            stageHead.style.background = '#B8B8B8';
          }
        }
      });
    }
    
    function getStageNameFromId(stageId) {
      const names = {
        'prune': 'Data pruning',
        'dataset_rep': 'Dataset reporting',
        'pheno_rep': 'Phenotype reporting',
        'contrast': 'Contrast selection',
        'ct': 'CT (convergence)',
        'ct_signif': 'CT signification (convergence)',
        'ct_disambig': 'CT disambiguation (convergence)',
        'ct_postproc': 'CT post-processing',
        'ora': 'ORA',
        'string': 'STRING',
        'ct_acc': 'CT accumulation (convergence)',
        'rer': 'RERconverge (RER)',
        'fade': 'FADE (selection)',
        'molerate': 'Molerate (RER)'
      };
      return names[stageId] || stageId;
    }
    
    function getStageColor(stageId) {
      const typeColors = {
        'prune': '#0EA5E9',
        'dataset_rep': '#7C3AED',
        'pheno_rep': '#7C3AED',
        'contrast': '#0EA5E9',
        'ct': '#F97316',
        'ct_signif': '#7C3AED',
        'ct_disambig': '#F97316',
        'ct_postproc': '#0EA5E9',
        'ora': '#7C3AED',
        'string': '#7C3AED',
        'ct_acc': '#F97316',
        'rer': '#F97316',
        'fade': '#F97316',
        'molerate': '#F97316'
      };
      return typeColors[stageId] || '#B8B8B8';
    }
  </script>
</body>
</html>
""".stripIndent()
    }

    // ── Dynamic directory scanning ─────────────────────────────────────────────

    static Map getCanonicalDirectories() {
        return [
            prune        : ['data_exploration/0.Data-pruning'],
            dataset_rep  : ['data_exploration'],
            pheno_rep    : ['data_exploration'],
            contrast     : ['data_exploration/2.CT', 'data_exploration/2.CT/1.Traitfiles'],
            ct           : ['caastools', 'discovery'],
            ct_signif    : ['signification', 'signification/gene_lists'],
            ct_disambig  : ['ct_disambiguation'],
            ct_postproc  : ['postproc', 'postproc/preprocessed'],
            ora          : ['ora'],
            string       : ['string'],
            ct_acc       : ['accumulation', 'accumulation/aggregation'],
            rer          : ['rerconverge'],
            fade         : ['selection/fade'],
            molerate     : ['selection/molerate']
        ]
    }

    static boolean checkStageCompletion(String baseDir, String stageId) {
        def canonicalDirs = getCanonicalDirectories()
        def requiredDirs = canonicalDirs[stageId]
        if (!requiredDirs) return false
        
        return requiredDirs.any { dirPath ->
            new File(baseDir, dirPath).exists()
        }
    }

    static Map scanWorkflowDirectory(String baseDir) {
        def canonicalDirs = getCanonicalDirectories()
        def results = [:]
        
        canonicalDirs.each { stageId, dirPaths ->
            results[stageId] = checkStageCompletion(baseDir, stageId)
        }
        
        return results
    }

    // ── Convenience: build ctx map from dynamic scanning or workflow context ───
    // Can be called with just a directory path for standalone operation,
    // or with params + workflow for backward compatibility.

    static Map buildCtx(baseDir, params = null, workflow = null) {
        def outdir = baseDir?.toString() ?: 
                    (params?.outdir ? params.outdir.toString() : "${workflow?.projectDir}/Out")
        
        // Dynamic scanning of workflow completion status
        def scanResults = scanWorkflowDirectory(outdir)
        
        [
            outdir        : outdir,
            launchDir     : workflow?.launchDir?.toString() ?: outdir,
            projectDir    : workflow?.projectDir?.toString() ?: outdir,
            profile       : workflow?.profile ?: 'standalone',
            runName       : workflow?.runName ?: 'dynamic-scan',
            sessionId     : workflow?.sessionId?.toString() ?: UUID.randomUUID().toString(),
            commandLine   : workflow?.commandLine ?: 'dynamic generation',
            prune         : scanResults.prune ?: false,
            datasetReport : scanResults.dataset_rep ?: false,
            phenotypeRep  : scanResults.pheno_rep ?: false,
            contrastSel   : scanResults.contrast ?: false,
            ct            : scanResults.ct ?: false,
            ctSignif      : scanResults.ct_signif ?: false,
            ctDisambig    : scanResults.ct_disambig ?: false,
            ctPostproc    : scanResults.ct_postproc ?: false,
            ora           : scanResults.ora ?: false,
            string        : scanResults.string ?: false,
            ctAccum       : scanResults.ct_acc ?: false,
            rer           : scanResults.rer ?: false,
            fade          : scanResults.fade ?: false,
            molerate      : scanResults.molerate ?: false
        ]
    }
}
