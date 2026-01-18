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
    tuple val(alignmentID), file(alignmentFile), path(discoveryFile, stageAs: 'discovery_*')
    path(resampledPath)  // Can be either directory or file
    
    output:
    tuple val(alignmentID), file("${alignmentID}.bootstraped.output"), optional: true

    script:
    def args = task.ext.args ?: ''
    // Allow discovery input from this run or a legacy single file / directory provided via params.discovery_out
    def discoveryCandidate = discoveryFile ?: (params.discovery_out != "none" ? file(params.discovery_out) : null)
    def discovery_arg = ''

    if (discoveryCandidate) {
        def candidatePath = discoveryCandidate instanceof java.nio.file.Path ? discoveryCandidate : file(discoveryCandidate)

        if (java.nio.file.Files.isDirectory(candidatePath)) {
            def nested = candidatePath.resolve("${alignmentID}.output")
            if (java.nio.file.Files.exists(nested)) {
                discovery_arg = "--discovery ${nested}"
            }
        } else if (java.nio.file.Files.exists(candidatePath)) {
            discovery_arg = "--discovery ${candidatePath}"
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