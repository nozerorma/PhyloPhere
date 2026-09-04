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

    publishDir = [
        path: { "${params.outdir}/bootstrap" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        enabled: params.publish_intermediates
    ]

    input:
    tuple val(alignmentID), path(alignmentFile), path(discoveryFile), path(resampledPath) // resampledPath can be either directory or file
    file caas_config

    output:
    tuple val(alignmentID), file("${alignmentID}.bootstraped.output"), emit: bootstrap_out
    tuple val(alignmentID), file("${alignmentID}.bootstrap.groups.output"), emit: bootstrap_groups, optional: true
    tuple val(alignmentID), file("${alignmentID}.bootstrap.discovery.output"), emit: bootstrap_perm_discovery, optional: true

    script:
    def args = "--patterns ${params.patterns} ${params.miss_pair ? '--miss_pair' : ''} ${params.caap_mode ? '--caap_mode' : ''}"
    // FOP (Gap A): resample over the fanned fop_labelings.tab (staged inside
    // resampledPath) and collapse recovery_boot to base-cycle units.
    def fop_on = (params.caas_perms_fop && params.multi_hypothesis)
    def fop_arg = fop_on ? "--fop" : ""
    def discovery_arg = discoveryFile.name != 'NO_FILE' ? "--discovery ${discoveryFile}" : ""
    // Belt-and-braces: --export_groups opens groups_handle, which forces boot.py
    // onto the scalar caasboot loop (use_vectorized requires `groups_handle is
    // None`). Under FOP that loop walks ~max_fop x cycles labelings per position,
    // so the flag turns a minutes-long vectorized run into a multi-hour one to
    // produce a ~90k-labeling groups dump nothing downstream consumes. Never emit
    // groups in FOP mode; the flag stays honored for plain-bootstrap debugging.
    def export_groups_arg = (!fop_on && params.export_groups != null && params.export_groups != "none") ? "--export_groups ${alignmentID}.bootstrap.groups.output" : ""
    def export_perm_discovery_arg = (params.export_perm_discovery != null && params.export_perm_discovery != "none") ? "--export_perm_discovery ${alignmentID}.bootstrap.discovery.output" : ""

    def pairArgs = """
if [ -d "${caas_config}" ]; then
    sample_file=\$(find -L ${caas_config} -type f -name '*.tab' | head -n 1)
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$sample_file" | sort -nu | wc -l | tr -d ' ')
else
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' ${caas_config} | sort -nu | wc -l | tr -d ' ')
fi
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
        # BLAS thread sizing: the vectorized CAAS kernel (modules/boot_vec.py) does
        # Level-3 matmuls; size OpenBLAS to the task allocation (NOT a global pin —
        # RERConverge shares this OpenBLAS and needs multithreaded BLAS). MKL/NUMEXPR
        # are unused here and pinned to 1 to avoid stray oversubscription.
        export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        ${pairArgs}
        /usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${caas_config} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${export_groups_arg} \\
            ${export_perm_discovery_arg} \\
            ${fop_arg} \\
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
        export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        ${pairArgs}
        $baseDir/subworkflows/CT/local/ct bootstrap \\
            -a ${alignmentFile} \\
            -t ${caas_config} \\
            -s ${resampledPath} \\
            -o ${alignmentID}.bootstraped.output \\
            --fmt ${params.ali_format} \\
            ${discovery_arg} \\
            ${export_groups_arg} \\
            ${export_perm_discovery_arg} \\
            ${fop_arg} \\
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
    label 'process_boot_batched'

    publishDir = [
        path: { "${params.outdir}/bootstrap" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
        enabled: params.publish_intermediates
    ]

    input:
    tuple val(batchID), val(batchSize), val(batchManifestText), path(alignmentFiles, stageAs: 'alignments/*'), path(discoveryFiles, stageAs: 'discovery/*'), path(resampledPath)
    file caas_config

    output:
    path("*.bootstraped.output"), emit: bootstrap_out
    path("*.bootstrap.groups.output"), emit: bootstrap_groups, optional: true
    path("*.bootstrap.discovery.output"), emit: bootstrap_perm_discovery, optional: true

    script:
    def args = "--patterns ${params.patterns} ${params.miss_pair ? '--miss_pair' : ''} ${params.caap_mode ? '--caap_mode' : ''}"
    // FOP (Gap A): base-cycle-collapsed bootstrap over the fanned fop_labelings.tab.
    def fopFlag = (params.caas_perms_fop && params.multi_hypothesis) ? '1' : '0'
    // Belt-and-braces: --export_groups forces the scalar caasboot loop (see the
    // non-batched process above). It is prohibitively expensive over the FOP
    // fan and produces a debug artifact nothing reads, so suppress it whenever
    // FOP is active — the flag still works for plain-bootstrap debugging.
    def exportGroupsFlag = (fopFlag == '0' && params.export_groups != null && params.export_groups != "none") ? '1' : '0'
    def ctBinary = (params.use_singularity || params.use_apptainer)
        ? "/usr/local/bin/_entrypoint.sh $baseDir/subworkflows/CT/local/ct"
        : "$baseDir/subworkflows/CT/local/ct"

    """
# Batched bootstrap fans out genes across an internal worker pool
# (--workers). Each worker runs the vectorized kernel, so pin BLAS to 1 thread
# per worker — parallelism comes from the pool, not from intra-matmul threading
# (otherwise workers x BLAS-threads oversubscribes the allocation).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cat > ${batchID}.manifest.tsv <<'EOF'
""" + batchManifestText + """EOF

cat > .ct_bootstrap_batch_args <<'EOF'
${args.replaceAll('\n', ' ')}
EOF

if [ -d "${caas_config}" ]; then
    sample_file=\$(find -L ${caas_config} -type f -name '*.tab' | head -n 1)
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' "\$sample_file" | sort -nu | wc -l | tr -d ' ')
else
    n_pairs=\$(awk '\$3~/^[0-9]+\$/{print \$3}' ${caas_config} | sort -nu | wc -l | tr -d ' ')
fi
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
    --workers ${task.cpus} \\
    --ali-format ${params.ali_format} \\
    --ct-bin ${ctBinary} \\
    --progress-log 0 \\
    --export-groups ${exportGroupsFlag} \\
    --export-perm-discovery ${params.export_perm_discovery != null && params.export_perm_discovery != "none" ? '1' : '0'} \\
    --fop ${fopFlag} \\
    --extra-args-file .ct_bootstrap_batch_args
"""
}
