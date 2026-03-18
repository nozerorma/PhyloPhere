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
    file caas_config

    output:
    tuple val(alignmentID), path("${alignmentID}.output"), emit: discovery_out, optional: true
    path("${alignmentID}.background.tsv"), emit: background_out, optional: true

    script:
    // Define extra discovery arguments from params.file
    def args = task.ext.args ?: ''

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh ct discovery \\
        -a ${alignmentFile} \\
        -t ${caas_config} \\
        -o ${alignmentID}.output \\
        --background_output ${alignmentID}.background.tsv \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')}
        """
    } else {
        """
        echo "Running locally"
        $baseDir/subworkflows/CT/local/ct discovery \\
        -a ${alignmentFile} \\
        -t ${caas_config} \\
        -o ${alignmentID}.output \\
        --background_output ${alignmentID}.background.tsv \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')}
        """
    }
}

process DISCOVERY_BATCHED {
    tag "$batchID (${batchSize} genes)"
    label 'process_discovery'

    input:
    tuple val(batchID), val(batchSize), val(batchManifestText), path(alignmentFiles, stageAs: 'alignments/*')
    file caas_config

    output:
    path("*.output"), emit: discovery_out, optional: true
    path("*.background.tsv"), emit: background_out, optional: true

    script:
    def args = task.ext.args ?: ''
    def runnerMode = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? '/usr/local/bin/_entrypoint.sh'
        : "$baseDir/subworkflows/CT/local/ct"

    """
    cat > ${batchID}.manifest.tsv <<'EOF'
    ${batchManifestText}
    EOF

    cat > .ct_discovery_batch_args <<'EOF'
    ${args.replaceAll('\n', ' ')}
    EOF

    bash $baseDir/subworkflows/CT/local/scripts/run_ct_discovery_batch.sh \\
        --batch-id ${batchID} \\
        --manifest ${batchID}.manifest.tsv \\
        --caas-config ${caas_config} \\
        --workers ${params.ct_discovery_batch_workers} \\
        --ali-format ${params.ali_format} \\
        --runner-mode ${runnerMode} \\
        --ct-bin ${ctBinary} \\
        --extra-args-file .ct_discovery_batch_args
    """.stripIndent()
}
