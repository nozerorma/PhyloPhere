#!/usr/bin/env nextflow
// rer_trees.nf — Pruned gene trees and RERconverge master tree object.
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_TREES: runs rer_master_tree.R, which prunes every gene tree to the species of the
 *  trait file (optionally renaming tips through the tax_id table) and reads the pruned
 *  trees with RERconverge::readTrees.
 *
 *  Consumes:  trait file, gene trees (multi-Newick file; unnamed trees are called gene1,
 *             gene2, ...), tax_id table (NO_FILE to skip the renaming)
 *  Produces:  <gene_trees>.pruned.txt (gene name, tab, Newick),
 *             <gene_trees>.masterTree.output (RDS of the readTrees object)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Master tree ────────────────────────────────────────────────────────────────

process RER_TREES {
    tag "$gene_trees_file"

    label 'process_reporting'
    // All gene trees are loaded into R at once: conf/resources.config gives RER_TREES
    // 32 GB x attempt, so maxRetries 3 steps through 32, 64 and 96 GB.
    maxRetries 3

    publishDir path: "${params.outdir}/rerconverge/rer_objects", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path my_traitfile
    path gene_trees_file
    path tax_id_file

    output:
    path("${gene_trees_file}.masterTree.output"), emit: master_tree
    path("${gene_trees_file}.pruned.txt"),        emit: pruned_trees

    script:
    def args = task.ext.args ?: ''
    def pruned_trees_out = "${gene_trees_file}.pruned.txt"
    def masterTrees_out = "${gene_trees_file}.masterTree.output"
    
    def tax_id_arg = (tax_id_file.name != 'NO_FILE') ? "${tax_id_file}" : ''

    if (params.use_singularity) {
        """
        echo "Using Singularity"
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/rer_master_tree.R' \\
        ${ gene_trees_file } \\
        ${ my_traitfile } \\
        ${ params.sp_colname } \\
        ${ pruned_trees_out } \\
        ${ masterTrees_out } \\
        ${ tax_id_arg } \\
        $args
        """
    } else {
        """
        echo "Running locally"
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/rer_master_tree.R' \\
        ${ gene_trees_file } \\
        ${ my_traitfile } \\
        ${ params.sp_colname } \\
        ${ pruned_trees_out } \\
        ${ masterTrees_out } \\
        ${ tax_id_arg } \\
        $args
        """
    }
}
