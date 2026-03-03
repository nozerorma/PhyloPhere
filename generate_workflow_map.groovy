#!/usr/bin/env groovy
/*
 * generate_workflow_map.groovy
 *
 * Standalone script for generating a dynamic PhyloPhere workflow map HTML
 * from any results directory. This script enables manual generation of the
 * workflow map without running the full Nextflow pipeline.
 *
 * The generated HTML is fully dynamic - it scans its parent directory for
 * canonical workflow directories and updates completion status accordingly.
 * It includes JavaScript refresh capability for real-time status updates.
 *
 * Usage:
 *   groovy generate_workflow_map.groovy [target_directory]
 *
 * If no target_directory is specified, uses current working directory.
 * The workflow_map.html file will be created in the target directory.
 *
 * Examples:
 *   # Generate map in current directory
 *   groovy generate_workflow_map.groovy
 *
 *   # Generate map in specific results directory
 *   groovy generate_workflow_map.groovy /path/to/results/
 *
 *   # Generate map in PhyloPhere output directory
 *   groovy generate_workflow_map.groovy ./Out/
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 */

@Grab('org.codehaus.groovy:groovy-all:2.4.15')

// Load the WorkflowMap library (assuming this script is in the same directory as lib/)
def scriptDir = new File(getClass().protectionDomain.codeSource.location.path).parent
def libDir = new File(scriptDir, 'lib')
if (!libDir.exists()) {
    println "ERROR: lib/ directory not found. Ensure this script is in the PhyloPhere root directory."
    System.exit(1)
}
// Add lib to classpath
this.class.classLoader.addURL(libDir.toURI().toURL())

// Load WorkflowMap class
def loader = new GroovyClassLoader()
loader.addURL(libDir.toURI().toURL())
def workflowMapFile = new File(libDir, 'WorkflowMap.groovy')
if (!workflowMapFile.exists()) {
    println "ERROR: WorkflowMap.groovy not found in lib/ directory."
    System.exit(1)
}
def workflowMapClass = loader.parseClass(workflowMapFile)

def targetDir = args.length > 0 ? args[0] : System.getProperty("user.dir")
def targetFile = new File(targetDir)

if (!targetFile.exists()) {
    println "ERROR: Target directory does not exist: ${targetDir}"
    System.exit(1)
}

if (!targetFile.isDirectory()) {
    println "ERROR: Target path is not a directory: ${targetDir}"
    System.exit(1)
}

def absoluteTargetDir = targetFile.canonicalPath

try {
    println "Generating workflow map for directory: ${absoluteTargetDir}"
    
    // Use dynamic scanning to build context
    def ctx = workflowMapClass.buildCtx(absoluteTargetDir)
    
    // Generate HTML
    def html = workflowMapClass.buildWorkflowMapHtml(ctx)
    
    // Write to target directory
    def outputFile = new File(targetFile, 'workflow_map.html')
    outputFile.text = html
    
    println "✓ Workflow map generated successfully:"
    println "  Location: ${outputFile.absolutePath}"
    println "  Mode: Dynamic directory scanning"
    
    // Show detected stages
    def scanResults = workflowMapClass.scanWorkflowDirectory(absoluteTargetDir)
    def completedStages = scanResults.findAll { it.value }
    if (completedStages) {
        println "  Detected completed stages:"
        completedStages.each { stageId, _ ->
            println "    - ${stageId}"
        }
    } else {
        println "  No completed stages detected"
    }
    
    println ""
    println "The HTML file includes:"
    println "  • Dynamic completion detection based on directory structure"
    println "  • JavaScript refresh capability for real-time updates"
    println "  • Full portability - works when moved to any location"
    
} catch (Exception e) {
    println "ERROR: Failed to generate workflow map: ${e.message}"
    e.printStackTrace()
    System.exit(1)
}