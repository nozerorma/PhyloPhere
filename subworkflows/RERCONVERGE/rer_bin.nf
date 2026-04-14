#!/usr/bin/env nextflow

/*
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: rer_bin.nf
#
*/

/*
 * ────────────────────────────────────────────────────────────────────────────
 * RER_BIN — binary RERconverge correlation
 *
 * Inputs
 * ──────
 *   trait_file      : RData with trait_vector (0/1 named numeric vector)
 *   rer_master_tree : geneTrees RDS from RER_TREES
 *   rer_matrix      : RER matrix RDS from RER_MATRIX
 *
 * Outputs
 * ───────
 *   fg_paths        : foreground paths RDS  (char2path equivalent)
 *   binary_output   : binary correlation results RDS
 * ────────────────────────────────────────────────────────────────────────────
 */

process RER_BIN {
    tag "$rer_matrix"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
    path trait_file
    path rer_master_tree
    path rer_matrix

    output:
    path "${params.traitname}.fg_paths.output",  emit: fg_paths
    path "${params.traitname}.binary.output",    emit: binary_output

    script:
    def perm_batches    = params.rer_perm_batches    ?: 0
    def perms_per_batch = params.rer_perms_per_batch ?: 100
    def min_pos         = params.rer_min_pos         ?: 2
    def binary_clade    = params.rer_binary_clade    ?: 'all'

    if (params.use_singularity || params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/binary_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.fg_paths.output \\
        ${rer_matrix} \\
        ${params.traitname}.binary.output \\
        ${params.rer_minsp} \\
        ${min_pos} \\
        ${params.winsorizeRER} \\
        ${binary_clade} \\
        ${perm_batches} \\
        ${perms_per_batch}
        """
    } else {
        """
        echo "Running locally"
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/binary_rer.R' \\
        ${trait_file} \\
        ${rer_master_tree} \\
        ${params.traitname}.fg_paths.output \\
        ${rer_matrix} \\
        ${params.traitname}.binary.output \\
        ${params.rer_minsp} \\
        ${min_pos} \\
        ${params.winsorizeRER} \\
        ${binary_clade} \\
        ${perm_batches} \\
        ${perms_per_batch}
        """
    }
}
