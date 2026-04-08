#!/usr/bin/env nextflow

/*
 * PHYLOPHERE: Site-PGLS workflow
 *
 * Runs phylogenetic GLS association tests on post-processed CAAS sites and
 * generates a compact HTML report.
 */

include { SITE_PGLS }   from "${baseDir}/subworkflows/PGLS/pgls_run.nf"
include { PGLS_REPORT } from "${baseDir}/subworkflows/PGLS/pgls_report.nf"

workflow PGLS {
    take:
        caas_input
        trait_input
        tree_input

    main:
        assert params.traitname : "PGLS requires --traitname"

        def caas_ch
        if (caas_input) {
            caas_ch = caas_input
        } else if (params.pgls_caas_input) {
            def f = file(params.pgls_caas_input)
            assert f.exists() : "PGLS: CAAS input not found: ${params.pgls_caas_input}"
            caas_ch = Channel.value(f)
        } else {
            caas_ch = Channel.empty()
        }

        def trait_ch = (trait_input ?: Channel.empty())
            .ifEmpty {
                def explicit = params.my_traits
                assert explicit : "PGLS requires --my_traits"
                def pruned = file("${params.outdir}/data_exploration/0.Data-pruning/pruned_trait_file.tsv")
                def resolved = pruned.exists() ? pruned : file(explicit)
                assert resolved.exists() : "PGLS: trait file not found: ${resolved}"
                resolved
            }

        def tree_ch = (tree_input ?: Channel.empty())
            .ifEmpty {
                def explicit = params.tree
                assert explicit : "PGLS requires --tree"
                def pruned = file("${params.outdir}/data_exploration/0.Data-pruning/pruned_tree_file.nwk")
                def resolved = pruned.exists() ? pruned : file(explicit)
                assert resolved.exists() : "PGLS: tree file not found: ${resolved}"
                resolved
            }

        def pgls_out = SITE_PGLS(caas_ch, trait_ch, tree_ch)
        PGLS_REPORT(
            pgls_out.pgls_tsv,
            pgls_out.site_diag_tsv,
            pgls_out.intval_tsv,
            pgls_out.extremes_tsv,
            pgls_out.excess_tsv
        )

    emit:
        pgls_tsv      = pgls_out.pgls_tsv
        site_diag_tsv = pgls_out.site_diag_tsv
        intval_tsv    = pgls_out.intval_tsv
        extremes_tsv  = pgls_out.extremes_tsv
        excess_tsv    = pgls_out.excess_tsv
}
