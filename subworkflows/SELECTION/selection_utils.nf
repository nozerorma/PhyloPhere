#!/usr/bin/env nextflow

/*
 * SELECTION / shared utility processes
 *
 * EXTRACT_EXTREME_SPECIES – derive top/bottom species lists from the contrast
 *                           traitfile (species<TAB>contrast_group<TAB>pair) or,
 *                           for back-compat, a trait_stats.csv with global_label
 */


include { selectionQuote } from './selection_contract'

process EXTRACT_EXTREME_SPECIES {
    cache 'deep'
    tag "extract_extreme_species"

    publishDir path: "${params.outdir}/selection/species_sets",
               mode: 'copy', overwrite: true

    input:
    tuple path(stats_file), val(input_hash)
    val code_hash

    output:
    path "top_species.txt", emit: top_species
    path "bottom_species.txt", emit: bottom_species

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local/src"
    if (params.use_singularity || params.use_apptainer) {
        """
        echo '${code_hash}:${input_hash}' > inputs.sha256
        /usr/local/bin/_entrypoint.sh python -B ${selectionQuote("${local_dir}/extract_extreme_species.py")} \
            --stats-file ${selectionQuote(stats_file)} \
            --out-top "top_species.txt" \
            --out-bottom "bottom_species.txt"
        """
    } else {
        """
        echo '${code_hash}:${input_hash}' > inputs.sha256
        python -B ${selectionQuote("${local_dir}/extract_extreme_species.py")} \
            --stats-file ${selectionQuote(stats_file)} \
            --out-top "top_species.txt" \
            --out-bottom "bottom_species.txt"
        """
    }
}
