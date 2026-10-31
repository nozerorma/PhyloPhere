#!/usr/bin/env nextflow
include { rerQuote; rerHash; rerSources } from './rer_contract'

process RER_TREES {
    tag "RER_TREES|${params.traitname}"
    label 'process_reporting'
    errorStrategy 'terminate'
    cache 'deep'
    publishDir path: "${params.outdir}/rerconverge/rer_objects", mode: 'copy'

    input:
    path my_traitfile
    path gene_trees_file
    path tax_id_file

    val source_hash
    output:
    path "${gene_trees_file}.masterTree.output", emit: master_tree
    path "${gene_trees_file}.pruned.txt", emit: pruned_trees
    path "${gene_trees_file}.pruned.txt.audit.tsv", emit: audit

    script:
    def input_hash = rerHash([my_traitfile, gene_trees_file, tax_id_file])
    def rs = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh Rscript' : 'Rscript'

    """
    # RER sources: ${source_hash}; inputs: ${input_hash}
    ${rs} ${rerQuote("${baseDir}/subworkflows/RERCONVERGE/local/rer_master_tree.R")} ${rerQuote(gene_trees_file)} ${rerQuote(my_traitfile)} ${rerQuote(params.sp_colname)} ${rerQuote("${gene_trees_file}.pruned.txt")} ${rerQuote("${gene_trees_file}.masterTree.output")} ${rerQuote(tax_id_file)} ${rerQuote(params.traitname)}
    """
}
