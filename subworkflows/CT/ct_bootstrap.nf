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
    //label 'big_mem'

    
    input:
    tuple val(alignmentID), file(alignmentFile)
    file(resampledFile)
    
    output:
    tuple val(alignmentID), file("${alignmentID}.bootstraped.output"), optional: true

    script:
    def args = task.ext.args ?: ''

    """    
    /usr/local/bin/_entrypoint.sh ct bootstrap \\
        -a ${alignmentFile} \\
        -t ${params.traitfile} \\
        -s ${resampledFile} \\
        -o ${alignmentID}.bootstraped.output \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')}
    """
}
