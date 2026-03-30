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
    path nw_tree,     stageAs: 'nw_tree.nwk'
    path caas_config, stageAs: 'caas_config.tab'
    path trait_val,   stageAs: 'traitvalues.tab'

    output:
    path("${nw_tree.baseName}.resampled.output/")

    script:
    def args = task.ext.args ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"

        if [ ! -f "${nw_tree}" ]; then
            echo "[ERROR] Missing tree input: ${nw_tree}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi
        if [ ! -f "${caas_config}" ]; then
            echo "[ERROR] Missing caas config input: ${caas_config}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi
        if [ ! -f "${trait_val}" ]; then
            echo "[ERROR] Missing trait values input: ${trait_val}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi

        mkdir -p ${nw_tree.baseName}.resampled.output
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/CT/local/permulations.R' \\
        "${nw_tree}" \\
        "${caas_config}" \\
        ${params.cycles} \\
        ${params.perm_strategy} \\
        "${trait_val}" \\
        ${nw_tree.baseName}.resampled.output \\
        ${params.chunk_size} \\
        ${params.include_b0}
        """
    } else {
        """
        echo "Running locally"

        if [ ! -f "${nw_tree}" ]; then
            echo "[ERROR] Missing tree input: ${nw_tree}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi
        if [ ! -f "${caas_config}" ]; then
            echo "[ERROR] Missing caas config input: ${caas_config}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi
        if [ ! -f "${trait_val}" ]; then
            echo "[ERROR] Missing trait values input: ${trait_val}" >&2
            echo "[DEBUG] workdir:" >&2
            pwd >&2
            ls -la >&2
            exit 1
        fi

        mkdir -p ${nw_tree.baseName}.resampled.output
        Rscript \\
        '$baseDir/subworkflows/CT/local/permulations.R' \\
        "${nw_tree}" \\
        "${caas_config}" \\
        ${params.cycles} \\
        ${params.perm_strategy} \\
        "${trait_val}" \\
        ${nw_tree.baseName}.resampled.output \\
        ${params.chunk_size} \\
        ${params.include_b0}
        """
    }
}
