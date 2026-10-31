#!/usr/bin/env nextflow

/*
 * CT observed workflow
 * Scores the observed labeling from a discovery.tab that already exists (--discovery_from): master CSV and meta_caas
 * tables. A live run gets them from the permulation core instead (CAAS_CORE, CAAS_CORE_OBSERVED).
 */

include { CAAS_OBSERVED; observedInputHash } from '../subworkflows/CT_DISAMBIGUATION/ct_observed'

workflow CT_OBSERVED {
    take:
        discovery_in
        trait_file_in
        tree_file_in
        hyp_pairs_in   // Channel<path> or null: contrast_hypotheses_pairs.tsv of this run's contrast selection
        taxid_map_in   // value channel: tax_id map of the species (curated by NAME_CURATION, or params.tax_id) or NO_FILE

    main:
        // The observed design and tree: those of the integrated run when given, --caas_config and --tree otherwise.
        def upstream_trait = (trait_file_in ?: Channel.empty())
        def upstream_tree = (tree_file_in ?: Channel.empty())

        def trait_file = upstream_trait.ifEmpty {
            def trait_file_param = params.caas_config
            if (!trait_file_param) {
                error "CT observed scoring requires a trait file from CT/contrast_selection or --caas_config"
            }
            assert file(trait_file_param).exists() : "Observed design not found: ${trait_file_param}"
            file(trait_file_param)
        }

        def tree_file = upstream_tree.ifEmpty {
            def tree_file_param = params.tree
            if (!tree_file_param) {
                error "CT observed scoring requires a tree file from CT/contrast_selection or --tree"
            }
            assert file(tree_file_param).isFile() : "Observed tree not found: ${tree_file_param}"
            file(tree_file_param)
        }

        // contrast_hypotheses_pairs.tsv: per-(hypothesis, domain) PSS weights for the per-side FOP pooling; absent ->
        // equal-weight node pooling. Resolution:
        //   1. --ct_disambig_hypotheses_pairs / --scoring_hypotheses_pairs;
        //   2. integrated run (hyp_pairs_in given): the file emitted by this run's contrast selection, or NO_HYP_PAIRS
        //      for a single-hypothesis run;
        //   3. standalone run: auto-discover in outdir, else NO_HYP_PAIRS.
        def hp_param = params.ct_disambig_hypotheses_pairs ?: params.scoring_hypotheses_pairs ?: ''
        def hyp_pairs_file
        if (hp_param) {
            assert file(hp_param).isFile() : "Requested hypothesis weights file not found: ${hp_param}"
            hyp_pairs_file = Channel.value(file(hp_param))
        } else if (hyp_pairs_in != null) {
            hyp_pairs_file = hyp_pairs_in.ifEmpty(file('NO_HYP_PAIRS')).toList().map { items ->
                assert items.size() == 1 : 'CT_OBSERVED requires one weights input'
                items[0]
            }
        } else {
            hyp_pairs_file = trait_file.map { design ->
                def adjacent = file("${design}/contrast_hypotheses_pairs.tsv")
                def legacy = file("${params.outdir}/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv")
                design.isDirectory() && adjacent.isFile() ? adjacent : (legacy.isFile() ? legacy : file('NO_HYP_PAIRS'))
            }
        }

        def alignment = params.alignment ? file(params.alignment) : file('NO_ALIGNMENT')
        assert !params.alignment || alignment.isDirectory() : "Alignment directory not found: ${params.alignment}"
        def ensembl = params.gene_ensembl_file ? file(params.gene_ensembl_file) : file('NO_ENSEMBL')
        assert !params.gene_ensembl_file || ensembl.isFile() : "Gene whitelist not found: ${params.gene_ensembl_file}"
        def observed = CAAS_OBSERVED(discovery_in, trait_file.map { d -> tuple(observedInputHash(d), d) }, tree_file, hyp_pairs_file,
                                    taxid_map_in ?: Channel.value(file('NO_FILE')),
                                    Channel.value(tuple(observedInputHash(alignment), alignment)), Channel.value(ensembl))

    emit:
        results_dir = observed.results_dir
        master_csv = observed.master_csv
        meta_caas_dir = observed.meta_caas_dir
        meta_caas = observed.meta_caas
        design = trait_file   // the design and tree this scoring resolved, for the steps that explain its results
        tree = tree_file
}
