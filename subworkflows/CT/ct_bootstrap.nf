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
    tuple val(alignmentID), path(alignmentFile), path(discoveryFile), path(resampledPath) // resampledPath can be either directory or file
    file caas_config
    
    output:
    tuple val(alignmentID), file("${alignmentID}.bootstraped.output"), emit: bootstrap_out
    tuple val(alignmentID), file("${alignmentID}.bootstrap.groups.output"), emit: bootstrap_groups, optional: true
    tuple val(alignmentID), file("${alignmentID}.bootstrap.discovery.output"), emit: bootstrap_perm_discovery, optional: true

    script:
    def args = task.ext.args ?: ''
    def discovery_arg = discoveryFile.name != 'NO_FILE' ? "--discovery ${discoveryFile}" : ""
    def progress_log_arg = params.progress_log != "none" ? "--progress_log ${alignmentID}.progress.log" : ""
    def export_groups_arg = (params.export_groups != null && params.export_groups != "none") ? "--export_groups ${alignmentID}.bootstrap.groups.output" : ""
    def export_perm_discovery_arg = (params.export_perm_discovery != null && params.export_perm_discovery != "none") ? "--export_perm_discovery ${alignmentID}.bootstrap.discovery.output" : ""

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
        /usr/local/bin/_entrypoint.sh ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${caas_config} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${progress_log_arg} \\
            ${export_groups_arg} \\
            ${export_perm_discovery_arg} \\
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
        $baseDir/subworkflows/CT/local/ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${caas_config} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${progress_log_arg} \\
            ${export_groups_arg} \\
            ${export_perm_discovery_arg} \\
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

process BOOTSTRAP_BATCHED {
    tag "$batchID (${batchSize} genes)"
    label 'process_boot'

    input:
    tuple val(batchID), val(batchSize), val(batchManifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(discoveryFiles, stageAs: 'discovery/*'), path(resampledPath)
    file caas_config

    output:
    path("*.bootstraped.output"), emit: bootstrap_out
    path("*.bootstrap.groups.output"), emit: bootstrap_groups, optional: true
    path("*.bootstrap.discovery.output"), emit: bootstrap_perm_discovery, optional: true

    script:
    def args = task.ext.args ?: ''
    def runnerMode = (params.use_singularity || params.use_apptainer) ? 'container' : 'local'
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? '/usr/local/bin/_entrypoint.sh'
        : "$baseDir/subworkflows/CT/local/ct"

    """
cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

cat > .ct_bootstrap_batch_args <<'EOF'
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
echo " --max_conserved \$_max_conserved --max_bg_gaps \$_max_bg_gaps --max_fg_gaps \$_max_fg_gaps --max_gaps \$_max_gaps --max_bg_miss \$_max_bg_miss --max_fg_miss \$_max_fg_miss --max_miss \$_max_miss" >> .ct_bootstrap_batch_args
echo "Resolved thresholds for \$n_pairs pairs: max_conserved=\$_max_conserved bg_gaps=\$_max_bg_gaps fg_gaps=\$_max_fg_gaps gaps=\$_max_gaps bg_miss=\$_max_bg_miss fg_miss=\$_max_fg_miss miss=\$_max_miss"

bash $baseDir/subworkflows/CT/local/scripts/run_ct_bootstrap_batch.sh \\
    --batch-id ${batchID} \\
    --manifest ${batchID}.manifest.tsv \\
    --caas-config ${caas_config} \\
    --resampled-path ${resampledPath} \\
    --workers ${params.ct_bootstrap_batch_workers} \\
    --ali-format ${params.ali_format} \\
    --runner-mode ${runnerMode} \\
    --ct-bin ${ctBinary} \\
    --progress-log ${params.progress_log != "none" ? '1' : '0'} \\
    --export-groups ${params.export_groups != null && params.export_groups != "none" ? '1' : '0'} \\
    --export-perm-discovery ${params.export_perm_discovery != null && params.export_perm_discovery != "none" ? '1' : '0'} \\
    --extra-args-file .ct_bootstrap_batch_args
"""
}
