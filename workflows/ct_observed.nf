#!/usr/bin/env nextflow

/*
 * CT observed workflow
 * Scores the observed labeling from a discovery.tab that already exists (--discovery_from): master CSV and meta_caas
 * tables. A live run gets them from the permulation core instead (CAAS_CORE, CAAS_CORE_OBSERVED).
 */

include { CAAS_OBSERVED } from '../subworkflows/CT_DISAMBIGUATION/ct_observed'

workflow CT_OBSERVED {
    take:
        discovery_in
        trait_file_in
        tree_file_in
        hyp_pairs_in   // Channel<path> or null: contrast_hypotheses_pairs.tsv of this run's contrast selection

    main:
        // The observed design and tree: those of the integrated run when given, --caas_config and --tree otherwise.
        def upstream_trait = (trait_file_in ?: Channel.empty())
        def upstream_tree = (tree_file_in ?: Channel.empty())

        def trait_file = upstream_trait.ifEmpty {
            def trait_file_param = params.caas_config
            if (!trait_file_param) {
                error "CT observed scoring requires a trait file from CT/contrast_selection or --caas_config"
            }
            file(trait_file_param)
        }

        def tree_file = upstream_tree.ifEmpty {
            def tree_file_param = params.tree
            if (!tree_file_param) {
                error "CT observed scoring requires a tree file from CT/contrast_selection or --tree"
            }
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
        if (hp_param && file(hp_param).exists()) {
            hyp_pairs_file = Channel.value(file(hp_param))
        } else if (hyp_pairs_in != null) {
            hyp_pairs_file = hyp_pairs_in.ifEmpty(file('NO_HYP_PAIRS')).first()
        } else {
            def auto = file("${params.outdir}/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv")
            hyp_pairs_file = Channel.value(auto.exists() ? auto : file('NO_HYP_PAIRS'))
        }

        def observed = CAAS_OBSERVED(discovery_in, trait_file, tree_file, hyp_pairs_file)

    emit:
        results_dir = observed.results_dir
        master_csv = observed.master_csv
        meta_caas = observed.meta_caas
        global_meta_caas = observed.global_meta_caas
}
