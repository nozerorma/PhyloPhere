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
# File: ct_discovery.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DISCOVERY module: This module is responsible for the discovery process based on input alignments.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


process DISCOVERY {
    tag "$alignmentID"
    label 'process_discovery'

    input:
    tuple val(alignmentID), file(alignmentFile)

    output:
    tuple val(alignmentID), file("${alignmentID}.output"), optional: true

    script:
    // Define extra discovery arguments from params.file
    def args = task.ext.args ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh ct discovery \\
        -a ${alignmentFile} \\
        -t ${params.traitfile} \\
        -o ${alignmentID}.output \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')}
        """
    } else {
        """
        echo "Running locally"
        $baseDir/subworkflows/CT/local/ct discovery \\
        -a ${alignmentFile} \\
        -t ${params.traitfile} \\
        -o ${alignmentID}.output \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')}
        """
    }
}
