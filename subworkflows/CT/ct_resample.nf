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
# File: ct_resample.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RESAMPLE module: This module is responsible for resampling based on different strategies.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


process RESAMPLE {
    tag "$nw_tree"

    // Uncomment the following lines to assign workload priority.
    label 'process_resample'


    input:
    path nw_tree
    path trait_val

    output:
    file("${nw_tree}.resampled.output")

    script:
    def args = task.ext.args ?: ''
    def strategyCommand = ""

    // Determine the strategy command based on the provided strategy
    if (params.strategy == "FGBG") {
        strategyCommand = "-f ${params.fgsize} -b ${params.bgsize} -m random"
    } else if (params.strategy == "TEMPLATE") {
        strategyCommand = "--bytemp ${params.template} -m random"
    } else if (params.strategy == "PHYLORESTRICTED") {
        strategyCommand = "--bytemp ${params.template} --limit_by_group ${params.bygroup} -m random"
    } else if (params.strategy == "BM") {
        strategyCommand = "--bytemp ${params.template} --traitvalues ${params.traitvalues} --mode bm --perm_strategy ${params.perm_strategy}"
    } else {
        exit 1, "Invalid strategy: ${params.strategy}"
    }

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/CT/local/permulations.R' \\
        ${nw_tree} \\
        ${params.traitfile} \\
        ${params.cycles} \\
        ${params.perm_strategy} \\
        ${trait_val} \\
        ${nw_tree}.resampled.output
        """
    } else {
        """
        echo "Running locally"
        Rscript \\
        '$baseDir/subworkflows/CT/local/permulations.R' \\
        ${nw_tree} \\
        ${params.traitfile} \\
        ${params.cycles} \\
        ${params.perm_strategy} \\
        ${trait_val} \\
        ${nw_tree}.resampled.output
        """
    }


    // Could probably be changed to so I can use the ct logic:
    // if (params.use_singularity = true) {
    // """
    //  /usr/local/bin/_entrypoint.sh ct resample \\
    //  -p ${nw_tree} \\
    //  -o ${nw_tree}.resampled.output \\
    //  ${strategyCommand} \\
    //  $args
    // """
    // } else
    // """
    // ct resample \\
    //  -p ${nw_tree} \\
    //  -o ${nw_tree}.resampled.output \\
    //  ${strategyCommand} \\
    //  $args
    // """
}
