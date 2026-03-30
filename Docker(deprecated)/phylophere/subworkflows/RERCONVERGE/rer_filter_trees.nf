#!/usr/bin/env nextflow

/*
 * FILTER_GENE_TREES
 * ─────────────────
 * Filters a RERConverge gene-trees file (tab-separated GENE_NAME\tNewick) to
 * retain only genes present in one or more gene-list files.  The two
 * directional gene lists from COLLECT_GENE_SETS (top + bottom) are unioned
 * before filtering so that both FG directions are covered in the single RER
 * analysis.
 *
 * Inputs
 * ──────
 *   gene_trees_file   : The full gene trees file (GENE\tNewick, one per line)
 *   gene_set_top      : Top-direction gene list (plain text, one gene per line)
 *   gene_set_bottom   : Bottom-direction gene list (plain text, one gene per line)
 *
 * Outputs
 * ───────
 *   filtered_trees    : Filtered gene trees file (same format as input)
 *
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 * File: subworkflows/RERCONVERGE/rer_filter_trees.nf
 */

process FILTER_GENE_TREES {
    tag "rer_filter_trees"
    label 'process_low'

    publishDir path: "${params.outdir}/rerconverge/gene_sets",
               mode: 'copy', overwrite: true,
               pattern: 'filtered_gene_trees.txt'

    input:
    path gene_trees_file
    path gene_set_top
    path gene_set_bottom

    output:
    path "filtered_gene_trees.txt", emit: filtered_trees

    script:
    def local_dir = "${baseDir}/subworkflows/RERCONVERGE/local"

    if (params.use_singularity || params.use_apptainer) {
        """
        # Merge top and bottom gene lists into one (union, deduplicated)
        cat "${gene_set_top}" "${gene_set_bottom}" \\
            | grep -v '^#' | grep -v '^[[:space:]]*\$' | sort -u > _combined_gene_list.txt

        /usr/local/bin/_entrypoint.sh python3 \\
            "${local_dir}/filter_gene_trees.py" \\
            "${gene_trees_file}" \\
            _combined_gene_list.txt \\
            filtered_gene_trees.txt
        """
    } else {
        """
        cat "${gene_set_top}" "${gene_set_bottom}" \\
            | grep -v '^#' | grep -v '^[[:space:]]*\$' | sort -u > _combined_gene_list.txt

        python3 \\
            "${local_dir}/filter_gene_trees.py" \\
            "${gene_trees_file}" \\
            _combined_gene_list.txt \\
            filtered_gene_trees.txt
        """
    }
}
