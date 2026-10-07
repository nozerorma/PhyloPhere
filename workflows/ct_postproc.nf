#!/usr/bin/env nextflow

/*
##
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: ct_postproc.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CT Post-Processing Workflow: Handles cluster filtering and characterization of
 *  CT discovery results. Provides parameter sweep (exploratory mode) and single
 *  parameter filtering (filter mode) with optional report generation.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local processes from subworkflows
include { CAAS_PREPARE_POSTPROC_INPUT; CT_FILTER; CT_FILTER_SUMMARY; CAAS_FILTER_GENES; CAAS_BACKGROUND_CLEANUP; clusterParameterGrid; clusterFileSuffix } from '../subworkflows/CT_POSTPROC/ctpp_clustfilter'
include { CT_POSTPROC_REPORT } from '../subworkflows/CT_POSTPROC/ctpp_characterization'
include { ASR_ROBUSTNESS } from './asr_robustness'

workflow CT_POSTPROC {
    take:
        disambiguation_input_channel      // Post-disambiguation master CSV (optional, can use --disambiguation_input instead)
        background_genes_channel     // Global background genes file from CT module (preferred)
        disambiguation_dir_channel   // Full ct_disambiguation/ directory for ASR robustness diagnostics (optional)
        hyp_pairs_channel            // contrast_hypotheses_pairs.tsv of this run's contrast selection (null: auto-discover in outdir)

    main:
        def filter_dir_ch = Channel.value("${params.outdir}/postproc")

        // Determine discovery file source: channel input or parameter
        def discovery_file_ch
        def discovery_file_obj
        
        // Check if using upstream outputs (integrated mode) or standalone mode
        if (disambiguation_input_channel) {
            log.info "📥 Using discovery file from upstream disambiguation output"
            discovery_file_ch = disambiguation_input_channel
            discovery_file_obj = null
        } else {
            assert params.disambiguation_input : "CT Post-Processing requires disambiguation master CSV from upstream workflow or --disambiguation_input parameter"
            discovery_file_obj = file(params.disambiguation_input)
            assert discovery_file_obj.exists() : "Error: disambiguation_input file not found: ${params.disambiguation_input}"
            assert discovery_file_obj.isFile() : "Error: disambiguation_input must be a file"
            discovery_file_ch = Channel.value(discovery_file_obj)
        }
        
        // ── ASR Robustness diagnostics (parallel, does NOT affect clustering path) ────────────
        // Resolve the ct_disambiguation/ directory: prefer the upstream channel, then
        // derive from disambiguation_input CSV path when running in standalone mode.
        if (params.asr_robustness) {
            def asr_dir_ch
            if (disambiguation_dir_channel) {
                asr_dir_ch = disambiguation_dir_channel
            } else if (params.disambiguation_input) {
                def csv_parent = file(params.disambiguation_input).parent
                asr_dir_ch = Channel.value(csv_parent)
            } else {
                asr_dir_ch = null
            }
            if (asr_dir_ch) {
                ASR_ROBUSTNESS(asr_dir_ch)
            } else {
                log.warn "[asr_robustness] No disambiguation directory available — ASR robustness skipped."
            }
        }

        // Normalize discovery/disambiguation schema and apply the current
        // precluster hard filters used by CT post-processing.
        // Contrast design for the species tally: this run's contrast selection, else the file in outdir, else the sentinel.
        def hyp_pairs_file
        if (hyp_pairs_channel != null) {
            hyp_pairs_file = hyp_pairs_channel.ifEmpty(file('NO_HYP_PAIRS')).first()
        } else {
            def auto_hp = file("${params.outdir}/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv")
            hyp_pairs_file = Channel.value(auto_hp.exists() ? auto_hp : file('NO_HYP_PAIRS'))
        }
        prepared_inputs = CAAS_PREPARE_POSTPROC_INPUT(discovery_file_ch, hyp_pairs_file)
        def prepared_discovery_ch = prepared_inputs.prepared_discovery

        log.info "📂 Post-processing input normalized from disambiguation master CSV"
        
        // Handle global background genes source for cleanup
        def global_background_genes
        if (params.ct_tool && background_genes_channel) {
            global_background_genes = background_genes_channel
            log.info "📥 Using global background genes from CT_DISCOVERY module"
        } else if (params.background_input) {
            def bg_path = file(params.background_input)
            assert bg_path.exists() : "Error: background_input file/directory not found: ${params.background_input}"

            if (bg_path.isDirectory()) {
                // Prefer background_genes-like files in standalone mode directories
                global_background_genes = Channel.fromPath("${params.background_input}/*background_genes*")
                log.info "📂 Loading global background genes from directory: ${params.background_input}"
            } else {
                // Single file
                global_background_genes = Channel.fromPath(params.background_input)
                log.info "📄 Loading global background genes file: ${params.background_input}"
            }
        } else {
            error "CT Post-Processing requires CT background_genes output or --background_input (global genes)"
        }
        
        // The parameter pair of the production filter. params.filter_minlen/filter_maxcaas may arrive as
        // plain Strings (JSON/CLI params), and CT_FILTER's script does `(maxcaas * 100).toInteger()`, where
        // `*` on a String repeats it ("0.7" * 100 is "0.70.7..."), so both are cast to numbers here.
        def filter_minlen_val = params.filter_minlen.toInteger()
        def filter_maxcaas_val = params.filter_maxcaas.toDouble()

        // Determine processing mode and create parameter combinations channel
        if (params.caas_postproc_mode == 'exploratory') {
            // Parameter sweep over minlen_values x maxcaas_values (plus the selected pair when it is outside the grid)
            def minlen_list = params.minlen_values.split(',').collect { it.trim().toInteger() }
            def maxcaas_list = params.maxcaas_values.split(',').collect { it.trim().toDouble() }
            def combos = clusterParameterGrid(minlen_list, maxcaas_list, filter_minlen_val, filter_maxcaas_val)

            param_combinations = Channel
                .fromList(combos)
                .combine(prepared_discovery_ch)
                .map { combo, disc_file -> tuple('exploratory', combo[0], combo[1], disc_file) }

            log.info "🔍 Exploratory mode: testing ${combos.size()} parameter combinations"
            log.info "   (Cluster filtering will respect caap_group boundaries)"

        } else if (params.caas_postproc_mode == 'filter') {
            // Single filter run: the selected minlen and maxcaas.
            param_combinations = Channel
                .of(tuple('filter', filter_minlen_val, filter_maxcaas_val))
                .combine(prepared_discovery_ch)
                .map { mode, minlen, maxcaas, disc_file ->
                    tuple(mode, minlen, maxcaas, disc_file)
                }
            
            log.info "🔧 Filter mode: running with minlen=${params.filter_minlen}, maxcaas=${params.filter_maxcaas}"
            log.info "   (Cluster filtering will respect caap_group boundaries)"
            
        } else {
            error "Invalid caas_postproc_mode: ${params.caas_postproc_mode}. Must be 'exploratory' or 'filter'"
        }
        
        // Run cluster filtering process
        filter_results = CT_FILTER(param_combinations)
        
        // Collect all results and generate consolidated summary
        filter_summary_results = CT_FILTER_SUMMARY(
            filter_results.filtered_files.collect()
        )
        
        // Run gene-level filtering (always -- see note below)
        def characterization_results = null

        // CAAS_FILTER_GENES always runs -- filter_caas_genes.py's own
        // mode='none' path (apply_gene_filter, mode='none') is a pure
        // passthrough of discovery_df, so this is the one code path that
        // correctly preserves the full discovery schema (asr_path_score,
        // core, ...) for every gene_filter_mode, 'none' included. An earlier
        // version special-cased 'none' to skip this process and reuse
        // CT_FILTER's own output instead, but that file only ever carries
        // Gene/Position/clustering_flag (see filter_caas_clusters-param.py's
        // docstring) -- not the annotated columns SCORING needs.
        assert params.gene_ensembl_file : "Error: --gene_ensembl_file is required for gene filtering"

        def gene_ensembl_file = file(params.gene_ensembl_file)
        assert gene_ensembl_file.exists() : "Error: gene_ensembl_file not found: ${params.gene_ensembl_file}"

        log.info "🧬 Running gene-level filtering (mode: ${params.gene_filter_mode})..."
        log.info "   (Gene-level statistics will be calculated per caap_group)"

        // The gene filter uses the cluster file of the selected pair (filter_minlen, filter_maxcaas),
        // in filter mode and in exploratory mode alike. The name follows CT_FILTER's output naming.
        def selected_cluster_suffix = clusterFileSuffix(filter_minlen_val, filter_maxcaas_val)
        def cluster_file = filter_results.filtered_files
            .collect()
            .map { files ->
                def hit = files.flatten().find { it.name.endsWith(selected_cluster_suffix) }
                assert hit : "Error: CT_FILTER produced no cluster file ending in ${selected_cluster_suffix}"
                hit
            }

        def gene_filter_results = CAAS_FILTER_GENES(
            prepared_discovery_ch,
            gene_ensembl_file,
            cluster_file
        )

        // Run background cleanup (always -- a no-op removed_genes_summary in
        // 'none' mode just yields an unchanged background)
        cleaned_backgrounds = CAAS_BACKGROUND_CLEANUP(
            global_background_genes,
            gene_filter_results.removed_genes
        )

        def filtered_discovery_ch = gene_filter_results.filtered_discovery
        def cleaned_background_main_ch = cleaned_backgrounds.cleaned_background_main

        log.info "Cleaned background files: ${params.outdir}/postproc/cleaned_backgrounds"
        
        // Characterization reports always run (gene_ensembl_file is already resolved and validated above).
        log.info "📊 CT characterization reports..."

        // The report reads the published cluster files of the filter mode's directory
        def filter_output_dir = params.caas_postproc_mode == 'exploratory' ?
            "${params.outdir}/postproc/filter_${params.caas_postproc_mode}" :
            "${params.outdir}/postproc/filter_selected"

        characterization_results = CT_POSTPROC_REPORT(
            prepared_discovery_ch,
            filter_summary_results.summary,
            filter_output_dir,
            gene_ensembl_file,
            gene_filter_results.gene_stats
        )

        log.info "Post-processing reports generated in: ${params.outdir}/postproc/reports"
    
    emit:
        filter_summary = filter_summary_results.summary  // filter_summary.tsv (rows=params, cols=groups)
        discarded_summary = filter_summary_results.discarded_summary  // discarded_summary.tsv (old format for compatibility)
        filter_dir = filter_dir_ch
        filtered_discovery = filtered_discovery_ch
        cleaned_background = cleaned_background_main_ch
}
