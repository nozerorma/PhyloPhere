#!/usr/bin/env nextflow

/*
 * EXTRACT_VEP_ARCHIVES
 * ────────────────────
 * Single-pass extraction of per-gene track (.html) and CDS FASTA files
 * from optional .tar.gz archives, keyed to the genes present in the CAAS file.
 *
 * When the inputs are plain directories the process creates a symlink so that
 * downstream processes (PROT2COORD, PROT2AA) always receive a real directory
 * regardless of whether the source was compressed or not.
 *
 * Two archive reads at most (one tar -tzf listing + one tar -xzf extraction)
 * rather than one read pair per gene.
 */

process EXTRACT_VEP_ARCHIVES {
    tag "extract_vep_archives"
    label 'process_long_compute'

    input:
    path caas_file
    path cds_src_in
    path track_src_in

    output:
    path "_cds_src",   emit: cds_dir
    path "_track_src", emit: track_dir

    script:
    """
    genes=\$(awk -F'\\t' '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (tolower(\$i) == "gene") { col=i; break }
            }
            next
        }
        col && \$col != "" { print \$col }
    ' "${caas_file}" | sort -u)
    gene_pat=\$(echo "\$genes" | paste -sd '|')

    if [[ "${cds_src_in}" == *.tar.gz || "${cds_src_in}" == *.tgz ]]; then
        mkdir -p _cds_src
        tar -tzf "${cds_src_in}" \
            | grep -E "(\${gene_pat})" | grep -E '\\.(fasta|fa|fna)\$' \
            > _cds_members.txt 2>/dev/null || true
        [[ -s _cds_members.txt ]] && tar -xzf "${cds_src_in}" -C _cds_src --files-from _cds_members.txt
    else
        ln -sf "\$(realpath "${cds_src_in}")" _cds_src
    fi

    if [[ "${track_src_in}" == *.tar.gz || "${track_src_in}" == *.tgz ]]; then
        mkdir -p _track_src
        tar -tzf "${track_src_in}" \
            | grep -E "(\${gene_pat})" | grep -E '\\.html\$' \
            > _track_members.txt 2>/dev/null || true
        [[ -s _track_members.txt ]] && tar -xzf "${track_src_in}" -C _track_src --files-from _track_members.txt
    else
        ln -sf "\$(realpath "${track_src_in}")" _track_src
    fi
    """
}
