#!/usr/bin/env nextflow
// ct_resample.nf — Harvest the pool of permulated FG/BG labelings used as the CAAS null.
// PhyloPhere | subworkflows/CT/

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
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RESAMPLE: runs permulations.R, which permulates the observed trait on the species
 *  tree and keeps params.caas_full_perms labelings that have the same number of
 *  independent pairs as the observed contrast (and, with params.multi_hypothesis, the
 *  same number of FOP hypotheses). The strategy is params.perm_strategy, the canonical pairs
 *  follow the PSS profile of the observed ones (params.perm_match_pss, tolerance
 *  params.perm_match_pss_tol; not for count traits) and the
 *  randomness is seeded with params.seed (1998 when unset).
 *
 *  Consumes:  species tree, observed traitfile (or a directory of traitfile_H*.tab,
 *             of which H1 is used), trait values
 *  Produces:  <tree>.resampled.output/ with resample_NNN.tab, permulation_manifest.tsv
 *             and, with multi_hypothesis, fop_labelings.tab and fop_pairs.tsv;
 *             copied to ${params.outdir}/resample when params.publish_intermediates
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Permulation harvest ──────────────────────────────────────────────────────

process RESAMPLE {
    tag "$nw_tree"

    label 'process_resample'

    publishDir path: { "${params.outdir}/resample" },
               mode: 'copy',
               saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
               enabled: params.publish_intermediates

    input:
    path nw_tree,     stageAs: 'nw_tree.nwk'
    path caas_config, stageAs: 'caas_config.tab'
    path trait_val,   stageAs: 'traitvalues/*'
    val species_column

    output:
    path("${nw_tree.baseName}.resampled.output/")

    script:
    def cycles = params.caas_full_perms == null ? 1000 : params.caas_full_perms
    if (!(cycles.toString() ==~ /[0-9]+/)) error "caas_full_perms must be a nonnegative integer."
    if (cycles.toString().toBigInteger() == 0) {
        def modulePath = groovy.json.JsonOutput.toJson("${baseDir}/subworkflows/CT/local")
        def py = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
        return """
        ${py} - "${caas_config}" "${nw_tree.baseName}.resampled.output" <<'EMPTY_HARVEST'
import json, sys
from pathlib import Path
sys.path.insert(0, ${modulePath})
from modules.design_io import read_design
read_design(sys.argv[1], ${params.multi_hypothesis == true ? 'True' : 'False'})
p = Path(sys.argv[2]); p.mkdir()
(p / 'permulation_provenance.json').write_text(json.dumps({'schema_version': 1, 'requested_cycles': 0, 'selected_cycles': 0, 'relaxed_reasons': []}) + '\\n')
EMPTY_HARVEST
        """
    }
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
        # With multi_hypothesis the input is a directory of traitfile_H*.tab; the observed pair
        # count comes from the canonical hypothesis H1; all active pairs are validated.
        if [ -d "${caas_config}" ]; then
            if [ -f "${caas_config}/traitfile_H1.tab" ]; then
                actual_caas_config="${caas_config}/traitfile_H1.tab"
            else
                echo "[ERROR] Canonical traitfile_H1.tab is missing" >&2
                exit 1
            fi
        else
            actual_caas_config="${caas_config}"
        fi
        if [ ! -f "\$actual_caas_config" ]; then
            echo "[ERROR] Missing caas config input: \$actual_caas_config" >&2
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
        '$baseDir/subworkflows/CT/local/scripts/permulations.R' \\
        "${nw_tree}" \\
        "\$actual_caas_config" \\
        ${params.caas_full_perms} \\
        ${params.perm_strategy} \\
        "${trait_val}" \\
        ${nw_tree.baseName}.resampled.output \\
        ${params.chunk_size} \\
        ${params.pss_top_pct != null ? params.pss_top_pct : 0.05} \\
        ${params.max_tries != null ? params.max_tries : 1000000} \\
        "${params.traitname ?: ''}" \\
        "${params.n_trait ?: ''}" \\
        "${params.c_trait ?: ''}" \\
        ${params.resample_use_n != null ? params.resample_use_n : true} \\
        "${params.trait_type ?: 'auto'}" \\
        ${params.multi_hypothesis ?: false} \\
        ${params.max_fop != null ? params.max_fop : 100} \\
        ${task.cpus} \\
        ${params.seed != null ? params.seed : 1998} \\
        ${params.perm_match_pss != null ? params.perm_match_pss : true} \\
        ${params.perm_match_pss_tol != null ? params.perm_match_pss_tol : 0.25} \\
        "${params.pair_center_gate ?: 'none'}" \
        ${params.ct_allow_relaxed_null == true} \
        "${species_column}"
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
        # With multi_hypothesis the input is a directory of traitfile_H*.tab; the observed pair
        # count comes from the canonical hypothesis H1; all active pairs are validated.
        if [ -d "${caas_config}" ]; then
            if [ -f "${caas_config}/traitfile_H1.tab" ]; then
                actual_caas_config="${caas_config}/traitfile_H1.tab"
            else
                echo "[ERROR] Canonical traitfile_H1.tab is missing" >&2
                exit 1
            fi
        else
            actual_caas_config="${caas_config}"
        fi
        if [ ! -f "\$actual_caas_config" ]; then
            echo "[ERROR] Missing caas config input: \$actual_caas_config" >&2
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
        '$baseDir/subworkflows/CT/local/scripts/permulations.R' \\
        "${nw_tree}" \\
        "\$actual_caas_config" \\
        ${params.caas_full_perms} \\
        ${params.perm_strategy} \\
        "${trait_val}" \\
        ${nw_tree.baseName}.resampled.output \\
        ${params.chunk_size} \\
        ${params.pss_top_pct != null ? params.pss_top_pct : 0.05} \\
        ${params.max_tries != null ? params.max_tries : 1000000} \\
        "${params.traitname ?: ''}" \\
        "${params.n_trait ?: ''}" \\
        "${params.c_trait ?: ''}" \\
        ${params.resample_use_n != null ? params.resample_use_n : true} \\
        "${params.trait_type ?: 'auto'}" \\
        ${params.multi_hypothesis ?: false} \\
        ${params.max_fop != null ? params.max_fop : 100} \\
        ${task.cpus} \\
        ${params.seed != null ? params.seed : 1998} \\
        ${params.perm_match_pss != null ? params.perm_match_pss : true} \\
        ${params.perm_match_pss_tol != null ? params.perm_match_pss_tol : 0.25} \\
        "${params.pair_center_gate ?: 'none'}" \
        ${params.ct_allow_relaxed_null == true} \
        "${species_column}"
        """
    }
}
