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

    def pairArgs = """
n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' ${caas_config} | sort -nu | wc -l | tr -d ' ')
_max_conserved=\$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}')
_max_bg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_bg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
echo "Resolved thresholds for \$n_pairs pairs: max_conserved=\$_max_conserved bg_gaps=\$_max_bg_gaps fg_gaps=\$_max_fg_gaps gaps=\$_max_gaps bg_miss=\$_max_bg_miss fg_miss=\$_max_fg_miss miss=\$_max_miss"
"""

    if (params.use_singularity | params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        ${pairArgs}
        /usr/local/bin/_entrypoint.sh ct discovery \\
        -a ${alignmentFile} \\
        -t ${caas_config} \\
        -o ${alignmentID}.output \\
        --background_output ${alignmentID}.background.tsv \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')} \\
        --max_conserved \$_max_conserved \\
        --max_bg_gaps \$_max_bg_gaps \\
        --max_fg_gaps \$_max_fg_gaps \\
        --max_gaps \$_max_gaps \\
        --max_bg_miss \$_max_bg_miss \\
        --max_fg_miss \$_max_fg_miss \\
        --max_miss \$_max_miss
        """
    } else {
        """
        echo "Running locally"
        ${pairArgs}
        $baseDir/subworkflows/CT/local/ct discovery \\
        -a ${alignmentFile} \\
        -t ${caas_config} \\
        -o ${alignmentID}.output \\
        --background_output ${alignmentID}.background.tsv \\
        --fmt ${params.ali_format} \\
        ${args.replaceAll('\n', ' ')} \\
        --max_conserved \$_max_conserved \\
        --max_bg_gaps \$_max_bg_gaps \\
        --max_fg_gaps \$_max_fg_gaps \\
        --max_gaps \$_max_gaps \\
        --max_bg_miss \$_max_bg_miss \\
        --max_fg_miss \$_max_fg_miss \\
        --max_miss \$_max_miss
        """
    }
}

process DISCOVERY_BATCHED {
    tag "$batchID (${batchSize} genes)"
    label 'process_discovery_batched'

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
""" + batchManifestText + """EOF

cat > .ct_discovery_batch_args <<'EOF'
${args.replaceAll('\n', ' ')}
EOF

n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' ${caas_config} | sort -nu | wc -l | tr -d ' ')
_max_conserved=\$(awk -v n="\$n_pairs" -v f="${params.min_divergent_fraction}" 'BEGIN{printf "%d", int(n*(1-f))}')
_max_bg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_gaps=\$(awk -v n="\$n_pairs" -v f="${params.max_gaps_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_bg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_bg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_fg_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_fg_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
_max_miss=\$(awk -v n="\$n_pairs" -v f="${params.max_miss_fraction}" 'BEGIN{printf "%d", int(n*f)}')
echo " --max_conserved \$_max_conserved --max_bg_gaps \$_max_bg_gaps --max_fg_gaps \$_max_fg_gaps --max_gaps \$_max_gaps --max_bg_miss \$_max_bg_miss --max_fg_miss \$_max_fg_miss --max_miss \$_max_miss" >> .ct_discovery_batch_args
echo "Resolved thresholds for \$n_pairs pairs: max_conserved=\$_max_conserved bg_gaps=\$_max_bg_gaps fg_gaps=\$_max_fg_gaps gaps=\$_max_gaps bg_miss=\$_max_bg_miss fg_miss=\$_max_fg_miss miss=\$_max_miss"

bash $baseDir/subworkflows/CT/local/scripts/run_ct_discovery_batch.sh \\
    --batch-id ${batchID} \\
    --manifest ${batchID}.manifest.tsv \\
    --caas-config ${caas_config} \\
    --workers ${params.ct_discovery_batch_workers} \\
    --ali-format ${params.ali_format} \\
    --runner-mode ${runnerMode} \\
    --ct-bin ${ctBinary} \\
    --extra-args-file .ct_discovery_batch_args
"""
}
