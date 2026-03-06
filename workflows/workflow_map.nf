/*
 * workflow_map.nf
 *
 * Dynamic PhyloPhere workflow-map HTML generator. Creates a self-contained HTML
 * artifact that scans its environment to determine workflow completion status.
 *
 * The generated HTML is fully portable and location-aware - it detects the 
 * presence of canonical workflow directories (caastools/, rerconverge/, etc.) 
 * to determine which stages have completed, regardless of where the HTML file
 * is moved within the results structure.
 *
 * **NEW**: Dynamic directory scanning replaces parameter-based detection.
 * The HTML includes JavaScript refresh capability to update status on demand.
 *
 * The **primary** map generation path is still the `workflow.onComplete` hook in
 * main.nf, which fires after ALL processes finish. Use this module when you need
 * the map inside the DAG as a publishDir-tracked output that depends on explicit
 * upstream channels.
 *
 * Usage example (downstream of a specific output):
 *
 *   include { WORKFLOW_MAP } from './workflows/workflow_map.nf'
 *   ...
 *   // Trigger after any results are ready:
 *   WORKFLOW_MAP(some_results.collect().map { true })
 *
 * The generated HTML can be moved anywhere within the results directory tree
 * and will dynamically adapt to its new location.
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

// ── Process ───────────────────────────────────────────────────────────────────

process GENERATE_WORKFLOW_MAP {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    val trigger

    output:
    path 'workflow_map.html'
    path 'workflow.html'

    exec:
    // Build context map using dynamic directory scanning instead of workflow parameters.
    // The baseDir is the output directory where the HTML will be placed.
    def outdir = params.outdir ? params.outdir.toString() : "${workflow.projectDir}/Out"
    def ctx = WorkflowMap.buildCtx(outdir, params, workflow)

    // Generate the HTML and write it to the process work directory so that
    // Nextflow can stage and publish the output file correctly.
    def html = WorkflowMap.buildWorkflowMapHtml(ctx)
    new File(task.workDir.toString(), 'workflow_map.html').text = html
    new File(task.workDir.toString(), 'workflow.html').text = html
}

// ── Workflow ──────────────────────────────────────────────────────────────────

workflow WORKFLOW_MAP {

    take:
    trigger   // val: any truthy value used as a synchronisation signal

    main:
    GENERATE_WORKFLOW_MAP(trigger)

    emit:
    html = GENERATE_WORKFLOW_MAP.out
}
