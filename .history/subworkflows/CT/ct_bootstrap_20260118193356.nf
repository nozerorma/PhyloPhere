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
    
    // Auto-detect discovery input type and construct appropriate argument
    def discovery_arg = ""
    
    if (discoveryFile) {
        // If discoveryFile is provided, detect if it's a directory or a file
        def discovery_path = "${discoveryFile}"
        def is_dir = new File(discovery_path).isDirectory()
        
        if (is_dir) {
            // If it's a directory, look for discovery files and use the one for this alignment
            def discovery_file_path = new File(discovery_path, "${alignmentID}.output")
            if (discovery_file_path.exists()) {
                discovery_arg = "--discovery ${discovery_file_path}"
            } else {
                println("Warning: No discovery file found for ${alignmentID} in directory ${discovery_path}")
            }
        } else {
            // It's a single file, use it directly
            discovery_arg = "--discovery ${discoveryFile}"
        }
    } else if (params.discovery_out != "none") {
        // Fallback to params.discovery_out with auto-detection
        def discovery_path = "${params.discovery_out}"
        def is_dir = new File(discovery_path).isDirectory()
        
        if (is_dir) {
            // It's a directory, look for the alignment-specific file
            def discovery_file_path = new File(discovery_path, "${alignmentID}.output")
            if (discovery_file_path.exists()) {
                discovery_arg = "--discovery ${discovery_file_path}"
            }
        } else {
            // It's a single file, use it directly
            discovery_arg = "--discovery ${discovery_path}"
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
