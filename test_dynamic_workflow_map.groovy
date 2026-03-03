#!/usr/bin/env groovy
/*
 * test_dynamic_workflow_map.groovy
 *
 * Test script for the dynamic workflow map functionality.
 * Tests directory scanning and completion detection logic.
 *
 * Usage:
 *   groovy test_dynamic_workflow_map.groovy
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

// Load the WorkflowMap library
def scriptDir = new File(getClass().protectionDomain.codeSource.location.path).parent
def libDir = new File(scriptDir, 'lib')
def loader = new GroovyClassLoader()
loader.addURL(libDir.toURI().toURL())
def workflowMapFile = new File(libDir, 'WorkflowMap.groovy')
def workflowMapClass = loader.parseClass(workflowMapFile)

// Test 1: Canonical directories mapping
println "=== Test 1: Canonical Directories Mapping ==="
def canonicalDirs = workflowMapClass.getCanonicalDirectories()
println "Canonical directories defined for ${canonicalDirs.size()} stages:"
canonicalDirs.each { stageId, dirPaths ->
    println "  ${stageId}: ${dirPaths.join(', ')}"
}
println ""

// Test 2: Stage completion detection with non-existent directory
println "=== Test 2: Stage Completion Detection (Non-existent) ==="
def testDir = "/tmp/nonexistent_phylophere_test"
def ctCompleted = workflowMapClass.checkStageCompletion(testDir, 'ct')
def rerCompleted = workflowMapClass.checkStageCompletion(testDir, 'rer')
println "CT completed in ${testDir}: ${ctCompleted}"
println "RER completed in ${testDir}: ${rerCompleted}"
println ""

// Test 3: Dynamic scanning of empty directory
println "=== Test 3: Dynamic Scanning (Empty Directory) ==="
def scanResults = workflowMapClass.scanWorkflowDirectory(testDir)
println "Scan results for ${testDir}:"
scanResults.each { stageId, completed ->
    println "  ${stageId}: ${completed ? 'COMPLETED' : 'not completed'}"
}
println ""

// Test 4: Create temporary directory structure and test detection
println "=== Test 4: Dynamic Detection with Mock Structure ==="
def tempDir = File.createTempDir("phylophere_test", "")
def caasToolsDir = new File(tempDir, "caastools")
def rerConvergeDir = new File(tempDir, "rerconverge")
def fadeDir = new File(tempDir, "selection/fade")

caasToolsDir.mkdirs()
rerConvergeDir.mkdirs()
fadeDir.mkdirs()

println "Created mock structure in: ${tempDir.absolutePath}"
println "  - caastools/ (for CT detection)"
println "  - rerconverge/ (for RER detection)"
println "  - selection/fade/ (for FADE detection)"

def mockScanResults = workflowMapClass.scanWorkflowDirectory(tempDir.absolutePath)
println "Mock scan results:"
mockScanResults.each { stageId, completed ->
    if (completed) {
        println "  ✓ ${stageId}: COMPLETED"
    } else {
        println "  ✗ ${stageId}: not completed"
    }
}

// Test 5: Context building
println ""
println "=== Test 5: Dynamic Context Building ==="
def ctx = workflowMapClass.buildCtx(tempDir.absolutePath)
println "Context built for: ${ctx.outdir}"
println "Profile: ${ctx.profile}"
println "Mode: ${ctx.runName}"
println "Detected completions:"
['ct', 'rer', 'fade'].each { stage ->
    def completed = ctx[stage] ?: ctx["${stage}Completed"] ?: false
    println "  ${stage}: ${completed ? 'YES' : 'NO'}"
}

// Test 6: HTML generation
println ""
println "=== Test 6: HTML Generation ==="
try {
    def html = workflowMapClass.buildWorkflowMapHtml(ctx)
    def outputFile = new File(tempDir, "test_workflow_map.html")
    outputFile.text = html
    println "✓ HTML generated successfully: ${outputFile.absolutePath}"
    println "  Size: ${html.length()} characters"
    println "  Contains JavaScript: ${html.contains('function refreshStatus()') ? 'YES' : 'NO'}"
    println "  Self-contained: ${html.contains('<script>') && html.contains('<style>') ? 'YES' : 'NO'}"
} catch (Exception e) {
    println "✗ HTML generation failed: ${e.message}"
}

// Cleanup
tempDir.deleteDir()
println ""
println "=== Tests Completed ==="
println "Temporary test directory cleaned up."