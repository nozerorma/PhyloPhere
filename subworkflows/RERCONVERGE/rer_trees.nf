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
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: rer_obj.R
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  DISCOVERY module: This module is responsible for the discovery process based on input alignments.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


process RER_TREES {
    tag "$gene_trees_file"

    // Uncomment the following lines to assign workload priority.
    label 'process_medium' // have to tell it that only if using cluster!!!!!!!


    input:
    path my_traitfile
    path gene_trees_file

    output:
    file("${gene_trees_file}.masterTree.output") into outputNameChannel
    file("${gene_trees_file}.pruned.txt") into prunedTreesChannel

    script:
    def args = task.ext.args ?: ''
    def pruned_trees_out = "${gene_trees_file}.pruned.txt"
    def masterTrees_out = "${gene_trees_file}.masterTree.output"

    """
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/rer_master_tree.R' \\
        ${ gene_trees_file } \\
        ${ my_traitfile } \\
        ${ params.sp_colname } \\
        ${ pruned_trees_out } \\
        ${ masterTrees_out } \\
        $args
    """
}
