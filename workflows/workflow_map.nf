/*
 * workflow_map.nf
 *
 * Optional tracked Nextflow process/workflow for generating the PhyloPhere
 * workflow-map HTML artifact within the DAG (useful when you need the map to
 * be chained after specific upstream outputs).
 *
 * The **primary** map generation path is the `workflow.onComplete` hook in
 * main.nf, which fires after ALL processes finish.  Use this module only when
 * you need the map inside the DAG as a publishDir-tracked output that depends
 * on explicit upstream channels.
 *
 * Usage example (downstream of a specific output):
 *
 *   include { WORKFLOW_MAP } from './workflows/workflow_map.nf'
 *   ...
 *   // Trigger only after accum_results is ready:
 *   WORKFLOW_MAP(accum_results.gene_lists.collect().map { true })
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

    exec:
    // Build context map using the shared WorkflowMap library (lib/WorkflowMap.groovy).
    def ctx = WorkflowMap.buildCtx(params, workflow)

    // Generate the HTML and write it to the process work directory so that
    // Nextflow can stage and publish the output file correctly.
    def html = WorkflowMap.buildWorkflowMapHtml(ctx)
    new File(task.workDir.toString(), 'workflow_map.html').text = html
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
