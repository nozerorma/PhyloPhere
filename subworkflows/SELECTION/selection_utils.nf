#!/usr/bin/env nextflow

/*
 * SELECTION / shared utility processes
 *
 * PHYLIP_TO_FASTA       – normalize a supported alignment file to FASTA
 * FILTER_FASTA_TO_TREE  – drop FASTA sequences whose taxa are absent from a tree
 * EXTRACT_EXTREME_SPECIES – derive top/bottom species lists from trait_stats.csv
 * ANNOTATE_TREE_FG        – label foreground leaf branches in a Newick tree for FADE
 * EXTRACT_FG_BRANCHES     – write a list of FG species names for MoleRate --branches args
 * COLLECT_GENE_SETS     – build TOP / BOTTOM gene lists from postproc outputs
 */


process EXTRACT_EXTREME_SPECIES {
    tag "extract_extreme_species"

    publishDir path: "${params.outdir}/selection/species_sets",
               mode: 'copy', overwrite: true

    input:
    path stats_file

    output:
    path "top_species.txt", emit: top_species
    path "bottom_species.txt", emit: bottom_species

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh python ${local_dir}/extract_extreme_species.py \
            --stats-file "${stats_file}" \
            --out-top "top_species.txt" \
            --out-bottom "bottom_species.txt"
        """
    } else {
        """
        python ${local_dir}/extract_extreme_species.py \
            --stats-file "${stats_file}" \
            --out-top "top_species.txt" \
            --out-bottom "bottom_species.txt"
        """
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// PHYLIP_TO_FASTA
// direction is carried through as a passthrough value so that downstream
// processes (ANNOTATE_TREE_FG, FADE_RUN, MOLERATE_RUN) receive a consistent
// (gene_id, direction, fasta) tuple without needing a separate join step.
// Existing FASTA inputs are copied through directly; PHYLIP-like inputs are
// converted to FASTA.
// ─────────────────────────────────────────────────────────────────────────────
process PHYLIP_TO_FASTA {
    tag "${gene_id}|${direction}"

    input:
    tuple val(gene_id), val(direction), path(alignment_file)

    output:
    tuple val(gene_id), val(direction), path("${gene_id}.fa"), emit: fasta

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    // Format detection and tar.gz extraction are handled at bash level so that
    // .tar.gz archives (where alignment_file.name would be the archive name, not
    // the member name) are handled correctly.
    def pyBin = (params.use_singularity || params.use_apptainer)
        ? '/usr/local/bin/_entrypoint.sh python'
        : 'python'
    """
    # Convert to FASTA (or copy directly if already FASTA).
    _lower=\$(echo "${alignment_file}" | tr '[:upper:]' '[:lower:]')
    if [[ "\$_lower" == *.fa || "\$_lower" == *.fasta ]]; then
        cp "${alignment_file}" "${gene_id}.fa"
    else
        ${pyBin} ${local_dir}/phylip_to_fasta.py "${alignment_file}" "${gene_id}.fa"
    fi
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// FILTER_FASTA_TO_TREE
// Restrict each per-gene FASTA to the phenotype-pruned tree taxa before any
// downstream HyPhy step sees the alignment. This prevents off-tree species
// from entering FADE / MoleRate simply because they are present in the raw
// alignment directory.
// ─────────────────────────────────────────────────────────────────────────────
process FILTER_FASTA_TO_TREE {
    tag "${gene_id}|${direction}"

    input:
    tuple val(gene_id), val(direction), path(fasta), path(tree)

    output:
    tuple val(gene_id), val(direction), path("${gene_id}.${direction}.tree_filtered.fa"), emit: fasta

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    def pyBin = (params.use_singularity || params.use_apptainer)
        ? '/usr/local/bin/_entrypoint.sh python'
        : 'python'
    """
    ${pyBin} ${local_dir}/filter_fasta_to_tree.py \\
        --tree "${tree}" \\
        --fasta "${fasta}" \\
        --output "${gene_id}.${direction}.tree_filtered.fa"
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// ANNOTATE_TREE_FG  (FADE: {Foreground} labels in Newick string)
// The FASTA alignment is required so that the tree can be pruned to only the
// taxa present in the alignment, preventing FADE from failing when the species
// tree has more tips than the per-gene alignment.
// --fasta_out produces a filtered FASTA containing only sequences present in
// the (possibly pruned) tree, which prevents the reverse mismatch — when the
// alignment has MORE sequences than tree tips (e.g. a taxon present in the
// alignment but absent from the species tree).
// ─────────────────────────────────────────────────────────────────────────────
process ANNOTATE_TREE_FG {
    tag "${gene_id}_${direction}"

    input:
    tuple val(gene_id), val(direction), path(fasta), path(fg_species_file), path(tree)

    output:
    tuple val(gene_id), val(direction), path("${gene_id}_${direction}_fg.nwk"), emit: annotated_tree, optional: true
    tuple val(gene_id), val(direction), path("${gene_id}_${direction}.fa"),     emit: filtered_fasta, optional: true

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh python ${local_dir}/annotate_tree_fg.py \
            --species-file "${fg_species_file}" \
            --tree      "${tree}" \
            --fasta     "${fasta}" \
            --fasta_out "${gene_id}_${direction}.fa" \
            --output    "${gene_id}_${direction}_fg.nwk"
        """
    } else {
        """
        python ${local_dir}/annotate_tree_fg.py \
            --species-file "${fg_species_file}" \
            --tree      "${tree}" \
            --fasta     "${fasta}" \
            --fasta_out "${gene_id}_${direction}.fa" \
            --output    "${gene_id}_${direction}_fg.nwk"
        """
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// EXTRACT_FG_BRANCHES  (MoleRate: plain list of FG taxon names)
// ─────────────────────────────────────────────────────────────────────────────
process EXTRACT_FG_BRANCHES {
    tag "${direction}"

    input:
    tuple val(direction), path(fg_species_file)

    output:
    tuple val(direction), path("fg_branches_${direction}.txt"), emit: fg_list

    script:
    def local_dir = "${baseDir}/subworkflows/SELECTION/local"
    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh python ${local_dir}/extract_fg_branches.py \
            --species-file "${fg_species_file}" \
            --output    "fg_branches_${direction}.txt"
        """
    } else {
        """
        python ${local_dir}/extract_fg_branches.py \
            --species-file "${fg_species_file}" \
            --output    "fg_branches_${direction}.txt"
        """
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// COLLECT_GENE_SETS
// Combines postproc TXT files into two gene lists.
// ─────────────────────────────────────────────────────────────────────────────
process COLLECT_GENE_SETS {
    tag "collect_gene_sets"

    publishDir path: "${params.outdir}/selection/gene_sets",
               mode: 'copy', overwrite: true

    input:
    path pp_top,     stageAs: 'pp_top.txt'
    path pp_bottom,  stageAs: 'pp_bottom.txt'

    output:
    path "gene_set_top.txt",    emit: gene_set_top
    path "gene_set_bottom.txt", emit: gene_set_bottom

    script:
    def local_dir  = "${baseDir}/subworkflows/SELECTION/local"
    // Build optional source arguments (only pass a file if it was actually staged)
    def pp_top_arg     = pp_top.name     != 'NO_FILE' ? "--postproc_top        pp_top.txt"     : ""
    def pp_bottom_arg  = pp_bottom.name  != 'NO_FILE' ? "--postproc_bottom     pp_bottom.txt"  : ""

    if (params.use_singularity || params.use_apptainer) {
        """
        /usr/local/bin/_entrypoint.sh python ${local_dir}/collect_gene_sets.py \
            ${pp_top_arg} \
            ${pp_bottom_arg} \
            --out_top    gene_set_top.txt \
            --out_bottom gene_set_bottom.txt
        """
    } else {
        """
        python ${local_dir}/collect_gene_sets.py \
            ${pp_top_arg} \
            ${pp_bottom_arg} \
            --out_top    gene_set_top.txt \
            --out_bottom gene_set_bottom.txt
        """
    }
}
