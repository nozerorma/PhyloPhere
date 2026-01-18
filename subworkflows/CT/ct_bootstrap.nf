#!/usr/bin/env nextflow

/*
#                          _              _
#                         | |            | |
#      ___ __ _  __ _ ___| |_ ___   ___ | |___
#    / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
#   | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#   \___\__,_|\__,_|___/\__\___/ \___/|_|___/
#
# A Convergent Amino Acid Substitution identification 
# and analysis toolbox
#
# Author:         Fabio Barteri (fabio.barteri@upf.edu)
# Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu),
#                 Xavier Farré (xfarrer@igtp.cat),
#                 David de Juan (david.juan@upf.edu),
#                 Miguel Ramon (miguel.ramon@upf.edu) - Nextflow Protocol Elaboration
#
# File: ct_bootstrap.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  BOOTSTRAP module: This module is responsible for bootstraping on different discovery groups.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

process BOOTSTRAP {
    tag "$alignmentID"
    label 'process_boot'
    
    input:
    tuple val(alignmentID), file(alignmentFile), path(discoveryFile)
    path(resampledPath)  // Can be either directory or file
    
    output:
    tuple val(alignmentID), file("${alignmentID}.bootstraped.output"), optional: true

    script:
    def args = task.ext.args ?: ''
    // Use discovery file from pipeline join first; fall back to params.discovery_out for standalone mode
    // discoveryFile will be an empty list [] when discovery doesn't run in the pipeline
    def discovery_arg = ''
    
    if (discoveryFile && discoveryFile.toString() != '[]') {
        // Discovery file from pipeline join - use directly
        discovery_arg = "--discovery ${discoveryFile}"
    } else if (params.discovery_out != "none") {
        // Standalone mode: discovery output from params
        def discoveryPath = file(params.discovery_out)
        if (discoveryPath.isDirectory()) {
            // Directory mode: look for alignment-specific output
            def nested = discoveryPath.resolve("${alignmentID}.output")
            if (nested.exists()) {
                discovery_arg = "--discovery ${nested}"
            }
        } else if (discoveryPath.exists()) {
            // Single file mode: use directly
            discovery_arg = "--discovery ${discoveryPath}"
        }
    }
    
    def progress_log_arg = params.progress_log != "none" ? "--progress_log ${alignmentID}.progress.log" : ""

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${params.traitfile} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${progress_log_arg} \\
            ${args.replaceAll('\n', ' ')}
        """
    } else {
        """    
        echo "Running locally"
        $baseDir/ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${params.traitfile} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${progress_log_arg} \\
            ${args.replaceAll('\n', ' ')}
        """
    }
}