/*
 * WorkflowMap.groovy
 *
 * Helper library for generating the final PhyloPhere workflow-map HTML artifact.
 * Placed in lib/ so it is automatically compiled and available to the main.nf
 * script's workflow.onComplete hook.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

import java.util.UUID
import java.util.Base64
import groovy.json.JsonOutput

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

    static String readTextIfExists(String path) {
        if (!path) return null
        def f = new File(path)
        return f.exists() ? f.getText('UTF-8') : null
    }

    static String imageDataUriIfExists(String path) {
        if (!path) return null
        def f = new File(path)
        if (!f.exists()) return null

        def lower = f.name.toLowerCase()
        def mime = lower.endsWith('.png') ? 'image/png'
                 : (lower.endsWith('.jpg') || lower.endsWith('.jpeg')) ? 'image/jpeg'
                 : 'application/octet-stream'

        def encoded = Base64.getEncoder().encodeToString(f.bytes)
        return "data:${mime};base64,${encoded}"
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

        def htmlCandidates = (st.htmlCandidates ?: []) as List
        def htmlRows = htmlCandidates ? htmlCandidates.collect { ht ->
            def file = new File(ht.toString())
            def exists = file.exists()
            def targetPath = ht.toString()
            if (!exists) {
                def parent = file.parentFile
                if (parent && parent.exists()) {
                    def base = file.name.replace(".html", "")
                    def matching = parent.listFiles().find { it.name.startsWith(base) && it.name.endsWith(".html") }
                    if (matching) {
                        exists = true
                        targetPath = matching.absolutePath
                    }
                }
            }
            def label  = exists ? 'available' : 'MISSING'
            def href   = relativeFromOutdir(outdir, targetPath)
            def text   = relativeFromOutdir(outdir, targetPath)
            return "<div class=\"link-row\">→ html: <a href=\"${href}\">${htmlEscape(text)}</a> <span class=\"status ${exists ? 'ok' : 'missing'}\">${label}</span></div>"
        }.join('\n') : ''

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
        ${htmlRows}
      </div>
    </div>
    """.stripIndent()
    }

    // ── Full HTML page ─────────────────────────────────────────────────────────

    static String buildWorkflowMapHtml(Map ctx) {
        def outdir     = ctx.outdir
        def projectDir = ctx.projectDir ?: outdir
        def logoDataUri = imageDataUriIfExists("${projectDir}/res/logo.png")

        def colors = [
            reporting : '#7C3AED',   // reports / Enrichment
            prepost   : '#0EA5E9',   // pre/post-processing / characterization
            processes : '#F97316'    // CT / RER / disambiguation / accumulation / selection tools
        ]

        def stages = [
            [ id: 'prune',       name: 'Data pruning',                    type: 'prepost',   ran: ctx.prune,
              filesDirs: ["${outdir}/data_exploration/0.Data-pruning"],
              htmlCandidates: ["${outdir}/html_reports/2.Phenotype_exploration_pruned.html"] ],

            [ id: 'dataset_rep', name: 'Dataset reporting',               type: 'reporting', ran: ctx.datasetReport,
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/html_reports/1.Dataset_exploration.html"] ],

            [ id: 'pheno_rep',   name: 'Phenotype reporting',             type: 'reporting', ran: ctx.phenotypeRep,
              filesDirs: ["${outdir}/data_exploration"],
              htmlCandidates: ["${outdir}/html_reports/2.Phenotype_exploration_complete.html"] ],

            // NOTE: directories.R (TRAIT_ANALYSIS) only ever creates 1.Traitfiles and
            // 2.Bootstrap_traitfiles under 2.CT — there is no 3.Tree subdirectory.
            [ id: 'contrast',    name: 'Contrast selection',              type: 'prepost',   ran: ctx.contrastSel,
              filesDirs: ["${outdir}/data_exploration/2.CT",
                          "${outdir}/data_exploration/2.CT/1.Traitfiles",
                          "${outdir}/data_exploration/2.CT/2.Bootstrap_traitfiles"],
              htmlCandidates: ["${outdir}/html_reports/3.CI-composition.html",
                                "${outdir}/html_reports/4.Independent_contrasts.html"] ],

            // discovery/resample/bootstrap are only published when
            // --publish_intermediates is set (default false, see conf/ct.config);
            // caastools/ (from CT_CONCAT) is unconditional and is what actually
            // drives the "ran" status below.
            [ id: 'ct',          name: 'CT (convergence)',                type: 'processes', ran: ctx.ct,
              filesDirs: ["${outdir}/caastools", "${outdir}/discovery",
                          "${outdir}/resample",  "${outdir}/bootstrap"],
              htmlCandidates: [] ],

            // NOTE: 7.CT_signification.Rmd only ever creates meta_caas/ (meta_dir <-
            // "meta_caas") — there is no gene_lists subdirectory under signification/.
            [ id: 'ct_signif',   name: 'CT signification (convergence)',  type: 'reporting', ran: ctx.ctSignif,
              filesDirs: ["${outdir}/signification",
                          "${outdir}/signification/meta_caas"],
              htmlCandidates: ["${outdir}/html_reports/7.CT_signification.html"] ],

            [ id: 'ct_disambig', name: 'CT disambiguation (convergence)', type: 'processes', ran: ctx.ctDisambig,
              filesDirs: ["${outdir}/ct_disambiguation"],
              htmlCandidates: [] ],

            [ id: 'asr_robustness', name: 'ASR Robustness',              type: 'reporting', ran: ctx.asrRobustness,
              filesDirs: ["${outdir}/asr_robustness"],
              htmlCandidates: ["${outdir}/html_reports/9.ASR_robustness.html"] ],

            [ id: 'ct_postproc', name: 'CT post-processing',              type: 'prepost',   ran: ctx.ctPostproc,
              filesDirs: ["${outdir}/postproc",
                          "${outdir}/postproc/preprocessed",
                          "${outdir}/postproc/gene_filtering",
                          "${outdir}/postproc/cleaned_backgrounds"],
              htmlCandidates: ["${outdir}/html_reports/8.Characterization_report.html"] ],

            [ id: 'ct_acc',      name: 'CT accumulation (convergence)',   type: 'processes', ran: ctx.ctAccum,
              filesDirs: ["${outdir}/accumulation",
                          "${outdir}/accumulation/aggregation",
                          "${outdir}/accumulation/top/randomization",
                          "${outdir}/accumulation/bottom/randomization",
                          "${outdir}/accumulation/all/randomization",
                          "${outdir}/accumulation/top/gene_lists",
                          "${outdir}/accumulation/bottom/gene_lists",
                          "${outdir}/accumulation/all/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/10.Accumulation_report.html"] ],

            [ id: 'vep',         name: 'VEP characterization',            type: 'prepost',   ran: ctx.vep,
              filesDirs: ["${outdir}/vep"],
              htmlCandidates: [] ],

            [ id: 'rer',         name: 'RERconverge (RER)',               type: 'processes', ran: ctx.rer,
              filesDirs: ["${outdir}/rerconverge",
                          "${outdir}/rerconverge/rer_results",
                          "${outdir}/rerconverge/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/5.RERconverge_report.html"] ],

            [ id: 'fade',        name: 'FADE (selection)',                type: 'processes', ran: ctx.fade,
              filesDirs: ["${outdir}/selection/fade",
                          "${outdir}/selection/fade/top",
                          "${outdir}/selection/fade/bottom",
                          "${outdir}/selection/fade/top/gene_lists",
                          "${outdir}/selection/fade/bottom/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/6.FADE_report_top.html",
                               "${outdir}/html_reports/6.FADE_report_bottom.html"] ],

            // NOTE: scoring_compute.R writes the 9 STRING-ready slices flat as
            // gene_lists/slice_*.tsv (no position/ or gene/ subdirectories, and no
            // overlap/ directory — that was a stale reference to a since-removed
            // layout).
            [ id: 'scoring',     name: 'CAAS Scoring',                    type: 'reporting', ran: ctx.scoring,
              filesDirs: ["${outdir}/scoring",
                          "${outdir}/scoring/gene_lists"],
              htmlCandidates: ["${outdir}/html_reports/11.Scoring_report.html"] ],

            // The main ENRICHMENT workflow (workflows/enrichment.nf) runs downstream of
            // --scoring when --enrichment is also set. It only ever runs per-module FCS
            // for CAAS ('fcs') and RER ('scoring/rer') — FADE and accumulation no longer
            // run their own FCS ranking (see subworkflows/ENRICHMENT/fcs.nf), so there is
            // no 12.FCS_fade/12.FCS_accumulation to look for. AMI is separate: the
            // centralized DOMINO-based run (--scoring_ami) lives under ENRICHMENT and
            // publishes to ami/*, while each of RER/FADE/accumulation can *also* run
            // its own AMI report directly from its own workflow file (--ami, a
            // distinct flag — see workflows/rerconverge.nf, fade.nf, ct_accumulation.nf),
            // publishing under that module's own output tree instead of ami/.
            [ id: 'fcs',         name: 'Functional enrichment (FCS)',     type: 'reporting', ran: ctx.fcs,
              filesDirs: ["${outdir}/fcs",
                          "${outdir}/scoring/rer"],
              htmlCandidates: ["${outdir}/html_reports/12.FCS_scoring.html",
                               "${outdir}/html_reports/12.FCS_rer.html"] ],

            [ id: 'ami',         name: 'AMI (DOMINO active modules)',     type: 'reporting', ran: ctx.ami,
              filesDirs: ["${outdir}/ami",
                          "${outdir}/ami/ami_results",
                          "${outdir}/ami/ami_summary",
                          "${outdir}/ami/ami_plots",
                          "${outdir}/ami/ami_networks",
                          // Per-module AMI (--ami, run from each module's own
                          // workflow file rather than centralized ENRICHMENT/DOMINO).
                          "${outdir}/rerconverge/ami",
                          "${outdir}/selection/fade/top/ami",
                          "${outdir}/selection/fade/bottom/ami",
                          "${outdir}/accumulation/top/ami",
                          "${outdir}/accumulation/bottom/ami",
                          "${outdir}/accumulation/all/ami"],
              htmlCandidates: ["${outdir}/html_reports/13.AMI_analysis.html",
                               "${outdir}/html_reports/13.AMI_rer.html",
                               "${outdir}/html_reports/13.AMI_fade.html",
                               "${outdir}/html_reports/13.AMI_accumulation.html"] ],

            [ id: 'posenrich', name: 'Position Enrichment',              type: 'reporting', ran: ctx.posenrich,
              filesDirs: ["${outdir}/posenrich",
                          "${outdir}/posenrich/gmts"],
              htmlCandidates: ["${outdir}/html_reports/14.Position_enrichment_report.html"] ],

            [ id: 'compare',     name: 'Cross-module comparison',         type: 'reporting', ran: ctx.compare,
              filesDirs: ["${outdir}/compare",
                          "${outdir}/compare/compare_results"],
              htmlCandidates: ["${outdir}/html_reports/15.Comparison_report.html"] ]

        ].collect { st -> st + [color: colors[st.type]] }

        def chainIds = ['prune','dataset_rep','pheno_rep','contrast','ct','ct_signif',
                        'ct_disambig','asr_robustness','ct_postproc','ct_acc','vep','rer','fade',
                        'scoring','fcs','ami','posenrich','compare']
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
            "${projectDir}/conf/ct_disambiguation.config",
            "${projectDir}/conf/ct_postproc.config",
            "${projectDir}/conf/ct_accumulation.config",
            "${projectDir}/conf/enrichment.config",
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
    <div class="subtitle">Complete chain from prune/reporting to scoring (always shown). Gray = not run, colored = run.</div>
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
    <span class="sw" style="background:${colors.processes};"></span> Analysis processes (CT, RER, FADE, disambiguation, accumulation)
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
      
      // Canonical directories for each stage — generated straight from
      // getCanonicalDirectories() (single source of truth shared with the
      // Groovy-side scan) rather than a hand-duplicated literal, so this
      // list can't silently drift out of sync with it again.
      const canonicalDirs = ${JsonOutput.toJson(getCanonicalDirectories())};

      // Map stage IDs to display stages
      const stageMapping = {
        'prune': 'prune',
        'dataset_rep': 'dataset_rep',
        'pheno_rep': 'pheno_rep',
        'contrast': 'contrast',
        'ct': 'ct',
        'ct_signif': 'ct_signif',
        'ct_disambig': 'ct_disambig',
        'asr_robustness': 'asr_robustness',
        'ct_postproc': 'ct_postproc',
        'ct_acc': 'ct_acc',
        'vep': 'vep',
        'rer': 'rer',
        'fade': 'fade',
        'scoring': 'scoring',
        'fcs': 'fcs',
        'ami': 'ami',
        'compare': 'compare',
        'posenrich': 'posenrich'
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
        'asr_robustness': 'ASR Robustness',
        'ct_postproc': 'CT post-processing',
        'ct_acc': 'CT accumulation (convergence)',
        'vep': 'VEP characterization',
        'rer': 'RERconverge (RER)',
        'fade': 'FADE (selection)',
        'scoring': 'CAAS Scoring',
        'fcs': 'Functional enrichment (FCS)',
        'ami': 'AMI (DOMINO active modules)',
        'compare': 'Cross-module comparison',
        'posenrich': 'Position Enrichment'
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
        'asr_robustness': '#7C3AED',
        'ct_postproc': '#0EA5E9',
        'ct_acc': '#F97316',
        'vep': '#0EA5E9',
        'rer': '#F97316',
        'fade': '#F97316',
        'scoring': '#7C3AED',
        'fcs': '#7C3AED',
        'ami': '#7C3AED',
        'compare': '#7C3AED',
        'posenrich': '#7C3AED'
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
            ct_signif    : ['signification', 'signification/meta_caas'],
            ct_disambig  : ['ct_disambiguation'],
            asr_robustness : ['asr_robustness'],
            ct_postproc  : ['postproc', 'postproc/preprocessed'],
            ct_acc       : ['accumulation', 'accumulation/aggregation'],
            vep          : ['vep'],
            rer          : ['rerconverge', 'rerconverge/rer_results'],
            fade         : ['selection/fade', 'selection/fade/top', 'selection/fade/bottom'],
            scoring      : ['scoring'],
            fcs          : ['fcs', 'scoring/rer'],
            ami          : ['ami', 'rerconverge/ami',
                            'selection/fade/top/ami', 'selection/fade/bottom/ami',
                            'accumulation/top/ami', 'accumulation/bottom/ami', 'accumulation/all/ami'],
            compare      : ['compare'],
            posenrich    : ['posenrich']
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
    // Explicit overloads are used instead of default parameters to avoid
    // MissingMethodException when Groovy's generated overloads are incomplete
    // due to classloader / compile-cache issues in Nextflow.

    static Map buildCtx(baseDir) {
        return buildCtx(baseDir, null, null)
    }

    static Map buildCtx(baseDir, paramsMap) {
        return buildCtx(baseDir, paramsMap, null)
    }

    static Map buildCtx(baseDir, paramsMap, workflowMeta) {
        def outdir = baseDir?.toString() ?:
                    (paramsMap?.outdir ? paramsMap.outdir.toString() : "${workflowMeta?.projectDir}/out")
        
        // Dynamic scanning of workflow completion status
        def scanResults = scanWorkflowDirectory(outdir)
        
        [
            outdir        : outdir,
            launchDir     : workflowMeta?.launchDir?.toString() ?: outdir,
            projectDir    : workflowMeta?.projectDir?.toString() ?: outdir,
            profile       : workflowMeta?.profile ?: 'standalone',
            runName       : workflowMeta?.runName ?: 'dynamic-scan',
            sessionId     : workflowMeta?.sessionId?.toString() ?: UUID.randomUUID().toString(),
            commandLine   : workflowMeta?.commandLine ?: 'dynamic generation',
            prune         : scanResults.prune ?: false,
            datasetReport : scanResults.dataset_rep ?: false,
            phenotypeRep  : scanResults.pheno_rep ?: false,
            contrastSel   : scanResults.contrast ?: false,
            ct            : scanResults.ct ?: false,
            ctSignif      : scanResults.ct_signif ?: false,
            ctDisambig    : scanResults.ct_disambig ?: false,
            asrRobustness : scanResults.asr_robustness ?: false,
            ctPostproc    : scanResults.ct_postproc ?: false,
            ctAccum       : scanResults.ct_acc ?: false,
            vep           : scanResults.vep ?: false,
            rer           : scanResults.rer ?: false,
            fade          : scanResults.fade ?: false,
            scoring       : scanResults.scoring ?: false,
            fcs           : scanResults.fcs ?: false,
            ami           : scanResults.ami ?: false,
            compare       : scanResults.compare ?: false,
            posenrich     : scanResults.posenrich ?: false
        ]
    }
}
