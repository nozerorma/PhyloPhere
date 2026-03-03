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
include { CAAS_PREPARE_POSTPROC_INPUT; CT_FILTER; CT_FILTER_SUMMARY; CAAS_FILTER_GENES; CAAS_BACKGROUND_CLEANUP } from "${baseDir}/subworkflows/CT_POSTPROC/ctpp_clustfilter"
include { CT_POSTPROC_REPORT } from "${baseDir}/subworkflows/CT_POSTPROC/ctpp_characterization"

workflow CT_POSTPROC {
    take:
        disambiguation_input_channel      // Post-disambiguation master CSV (optional, can use --disambiguation_input instead)
        background_files_channel     // Raw background files from CT module (optional)
        background_genes_channel     // Global background genes file from CT module (preferred)
        bootstrap_input_channel      // kept for API compatibility (unused)
    
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
        
        // Normalize discovery/disambiguation schema and remove ambiguous/no_change before clustering
        prepared_inputs = CAAS_PREPARE_POSTPROC_INPUT(discovery_file_ch)
        def prepared_discovery_ch = prepared_inputs.prepared_discovery
        def precluster_removed_ch = prepared_inputs.removed_patterns

        log.info "📂 Post-processing input normalized from disambiguation master CSV"
        log.info "⛔ Rows with Pattern in {ambiguous, no_change} are removed before clustering and tracked"
        
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
        
        // Determine processing mode and create parameter combinations channel
        if (params.caas_postproc_mode == 'exploratory') {
            // Parameter sweep: generate cartesian product of all parameter combinations
            def minlen_list = params.minlen_values.split(',').collect { it.trim().toInteger() }
            def maxcaas_list = params.maxcaas_values.split(',').collect { it.trim().toDouble() }
            
            // Create channel with all parameter combinations
            // Combine parameters with discovery file channel
            param_combinations = Channel
                .from(minlen_list)
                .combine(Channel.from(maxcaas_list))
                .combine(prepared_discovery_ch)
                .map { minlen, maxcaas, disc_file -> 
                    tuple('exploratory', minlen, maxcaas, disc_file)
                }
            
            log.info "🔍 Exploratory mode: testing ${minlen_list.size()} × ${maxcaas_list.size()} = ${minlen_list.size() * maxcaas_list.size()} parameter combinations"
            log.info "   (Cluster filtering will respect CAAP_Group boundaries)"
            
        } else if (params.caas_postproc_mode == 'filter') {
            // Single filter run: use provided minlen and maxcaas
            // Combine with discovery file channel
            param_combinations = Channel
                .of(tuple('filter', params.filter_minlen, params.filter_maxcaas))
                .combine(prepared_discovery_ch)
                .map { mode, minlen, maxcaas, disc_file ->
                    tuple(mode, minlen, maxcaas, disc_file)
                }
            
            log.info "🔧 Filter mode: running with minlen=${params.filter_minlen}, maxcaas=${params.filter_maxcaas}"
            log.info "   (Cluster filtering will respect CAAP_Group boundaries)"
            
        } else {
            error "Invalid caas_postproc_mode: ${params.caas_postproc_mode}. Must be 'exploratory' or 'filter'"
        }
        
        // Run cluster filtering process
        filter_results = CT_FILTER(param_combinations)
        
        // Collect all results and generate consolidated summary
        filter_summary_results = CT_FILTER_SUMMARY(
            filter_results.filtered_files.collect()
        )
        
        // Run gene-level filtering if enabled
        def gene_filter_results = null
        def characterization_results = null
        def filtered_discovery_ch = Channel.empty()
        def cleaned_background_main_ch = Channel.empty()
        def ora_gene_lists_files_ch = Channel.empty()
        
        if (params.gene_filter_mode != 'none') {
            assert params.gene_ensembl_file : "Error: --gene_ensembl_file is required for gene filtering"
            
            def gene_ensembl_file = file(params.gene_ensembl_file)
            assert gene_ensembl_file.exists() : "Error: gene_ensembl_file not found: ${params.gene_ensembl_file}"
            
            log.info "🧬 Running gene-level filtering (mode: ${params.gene_filter_mode})..."
            log.info "   (Gene-level statistics will be calculated per CAAP_Group)"
            
            // Select appropriate cluster file based on mode
            // Use first cluster file from results without `.first()` to avoid
            // value-channel operator warnings.
            def cluster_file = filter_results.filtered_files
                .collect()
                .map { files ->
                    assert files && files.size() > 0 : "Error: CT_FILTER produced no cluster files"
                    files[0]
                }
            
            gene_filter_results = CAAS_FILTER_GENES(
                prepared_discovery_ch,
                gene_ensembl_file,
                cluster_file
            )
            
            // Run background cleanup (always when gene filtering is enabled)            
            cleaned_backgrounds = CAAS_BACKGROUND_CLEANUP(
                global_background_genes,
                gene_filter_results.removed_genes
            )

            filtered_discovery_ch = gene_filter_results.filtered_discovery
            cleaned_background_main_ch = cleaned_backgrounds.cleaned_background_main
            
            log.info "Cleaned background files: ${params.outdir}/postproc/cleaned_backgrounds"
        } else {
            // No gene filtering, no cleaned backgrounds
            gene_filter_results = null
            filtered_discovery_ch = Channel.empty()
            cleaned_background_main_ch = Channel.empty()
        }
        
        // Run characterization if reports are enabled
        if (params.generate_reports) {
            assert params.gene_ensembl_file : "Error: --gene_ensembl_file is required when generate_reports is true"
            
            def gene_ensembl_file = file(params.gene_ensembl_file)
            assert gene_ensembl_file.exists() : "Error: gene_ensembl_file not found: ${params.gene_ensembl_file}"
            
            log.info "📊 CT characterization reports..."
            
            // Pass the filter_ch output directory path instead of individual files
            def filter_output_dir = params.caas_postproc_mode == 'exploratory' ? 
                "${params.outdir}/postproc/filter_${params.caas_postproc_mode}" :
                "${params.outdir}/postproc/filter_selected"
            
            characterization_results = CT_POSTPROC_REPORT(
                prepared_discovery_ch,
                filter_summary_results.summary,
                filter_output_dir,
                gene_ensembl_file,
                gene_filter_results ? gene_filter_results.gene_stats : Channel.empty()
            )

            // Only pipe the specific txt exports needed by ORA (us_gs_relations/exports/txt/*.txt)
            // NOTE: disambiguation_characterization can emit a collection of files; normalize to
            // single-file items first to avoid calling String methods on ArrayList values.
            ora_gene_lists_files_ch = characterization_results.disambiguation_characterization
                .flatMap { item ->
                    item instanceof Collection ? item : [item]
                }
                .filter { f ->
                    def p = f.toString()
                    p.endsWith('.txt') && p.contains('us_gs_relations/exports/txt')
                }
            
            log.info "Post-processing reports generated in: ${params.outdir}/postproc/reports"
        }
    
    emit:
        filter_summary = filter_summary_results.summary  // filter_summary.tsv (rows=params, cols=groups)
        discarded_summary = filter_summary_results.discarded_summary  // discarded_summary.tsv (old format for compatibility)
        filter_dir = filter_dir_ch
        precluster_removed_patterns = precluster_removed_ch
        filtered_discovery = filtered_discovery_ch
        cleaned_background = cleaned_background_main_ch
        ora_gene_lists_files = ora_gene_lists_files_ch
}
