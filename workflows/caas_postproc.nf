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
# File: caas_postproc.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CAAS Post-Processing Workflow: Handles cluster filtering and characterization of
 *  CAAS discovery results. Provides parameter sweep (exploratory mode) and single
 *  parameter filtering (filter mode) with optional report generation.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local processes from subworkflows
include { CAAS_FILTER; CAAS_FILTER_SUMMARY; CAAS_FILTER_GENES; CAAS_BACKGROUND_CLEANUP } from "${baseDir}/subworkflows/CAAS_POSTPROC/ctpp_clustfilter"
include { CAAS_POSTPROC_REPORT } from "${baseDir}/subworkflows/CAAS_POSTPROC/ctpp_characterization"

workflow CAAS_POSTPROC {
    take:
        discovery_input_channel      // Discovery file from CT module (optional, can use --discovery_input instead)
        background_files_channel     // Background files from CT module (optional, can use --background_input instead)
        background_genes_channel     // Background genes file from CT module (optional)
    
    main:
        // Determine discovery file source: channel input or parameter
        def discovery_file_ch
        def discovery_file_obj
        
        // Check if using CT outputs (when ct_tool is enabled) or standalone mode
        if (params.ct_tool && discovery_input_channel) {
            log.info "📥 Using discovery file from CT_DISCOVERY module"
            discovery_file_ch = discovery_input_channel
            // For CAAP detection, we'll check later when the file is available
            discovery_file_obj = null
        } else {
            assert params.discovery_input : "CAAS Post-Processing requires either CT discovery output or --discovery_input parameter"
            discovery_file_obj = file(params.discovery_input)
            assert discovery_file_obj.exists() : "Error: discovery_input file not found: ${params.discovery_input}"
            assert discovery_file_obj.isFile() : "Error: discovery_input must be a file"
            discovery_file_ch = Channel.value(discovery_file_obj)
        }
        
        // Detect CAAP mode by checking for CAAP_Group column in discovery file
        // Only check if we have a file object (standalone mode)
        def is_caap_mode = false
        def caap_groups = []
        
        if (discovery_file_obj) {
            // Read first few lines to check for CAAP_Group column
            discovery_file_obj.withReader { reader ->
                def header = reader.readLine()
                if (header && header.contains("CAAP_Group")) {
                    is_caap_mode = true
                    log.info "📂 CAAP mode detected: discovery.tab contains CAAP_Group column"
                } else {
                    log.info "📄 CAAS mode detected: standard discovery.tab format"
                }
            }
        } else {
            log.info "📄 CAAS mode: will detect CAAP format during processing"
        }
        
        // Handle background files: use channel input if provided, otherwise use parameter
        def background_files
        if (params.ct_tool && background_files_channel) {
            background_files = background_files_channel
            log.info "📥 Using background files from CT_DISCOVERY module"
        } else if (params.background_input) {
            def bg_path = file(params.background_input)
            assert bg_path.exists() : "Error: background_input file/directory not found: ${params.background_input}"
            
            if (bg_path.isDirectory()) {
                // Directory: collect all background files
                background_files = Channel.fromPath("${params.background_input}/*.background.txt")
                log.info "📂 Loading background files from directory: ${params.background_input}"
            } else {
                // Single file
                background_files = Channel.fromPath(params.background_input)
                log.info "📄 Loading background file: ${params.background_input}"
            }
        } else {
            error "CAAS Post-Processing requires either CT background output or --background_input parameter"
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
                .combine(discovery_file_ch)
                .map { minlen, maxcaas, disc_file -> 
                    tuple('exploratory', minlen, maxcaas, disc_file)
                }
            
            log.info "🔍 Exploratory mode: testing ${minlen_list.size()} × ${maxcaas_list.size()} = ${minlen_list.size() * maxcaas_list.size()} parameter combinations"
            if (is_caap_mode) {
                log.info "   (Cluster filtering will respect CAAP_Group boundaries)"
            }
            
        } else if (params.caas_postproc_mode == 'filter') {
            // Single filter run: use provided minlen and maxcaas
            // Combine with discovery file channel
            param_combinations = Channel
                .of(tuple('filter', params.filter_minlen, params.filter_maxcaas))
                .combine(discovery_file_ch)
                .map { mode, minlen, maxcaas, disc_file ->
                    tuple(mode, minlen, maxcaas, disc_file)
                }
            
            log.info "🔧 Filter mode: running with minlen=${params.filter_minlen}, maxcaas=${params.filter_maxcaas}"
            if (is_caap_mode) {
                log.info "   (Cluster filtering will respect CAAP_Group boundaries)"
            }
            
        } else {
            error "Invalid caas_postproc_mode: ${params.caas_postproc_mode}. Must be 'exploratory' or 'filter'"
        }
        
        // Run cluster filtering process
        filter_results = CAAS_FILTER(param_combinations)
        
        // Collect all results and generate consolidated summary
        filter_summary = CAAS_FILTER_SUMMARY(
            filter_results.filtered_files.collect()
        )
        
        // Run gene-level filtering if enabled
        def gene_filter_results = null
        def cleaned_backgrounds = Channel.empty()
        
        if (params.gene_filter_mode != 'none') {
            assert params.gene_ensembl_file : "Error: --gene_ensembl_file is required for gene filtering"
            
            def gene_ensembl_file = file(params.gene_ensembl_file)
            assert gene_ensembl_file.exists() : "Error: gene_ensembl_file not found: ${params.gene_ensembl_file}"
            
            log.info "🧬 Running gene-level filtering (mode: ${params.gene_filter_mode})..."
            if (is_caap_mode) {
                log.info "   (Gene-level statistics will be calculated per CAAP_Group)"
            }
            
            // Select appropriate cluster file based on mode
            // Use first cluster file from results
            def cluster_file = filter_results.filtered_files.map { file -> file }.first()
            
            gene_filter_results = CAAS_FILTER_GENES(
                discovery_file_ch,
                gene_ensembl_file,
                cluster_file
            )
            
            // Run background cleanup (always when gene filtering is enabled)
            log.info "🧹 Cleaning background gene lists..."
            
            cleaned_backgrounds = CAAS_BACKGROUND_CLEANUP(
                background_files,
                gene_filter_results.removed_genes
            )
            
            log.info "Cleaned background files: ${params.outdir}/postproc/cleaned_backgrounds"
        } else {
            // No gene filtering, no cleaned backgrounds
            gene_filter_results = null
            cleaned_backgrounds = Channel.empty()
        }
        
        // Run characterization if reports are enabled
        if (params.generate_reports) {
            assert params.gene_ensembl_file : "Error: --gene_ensembl_file is required when generate_reports is true"
            
            def gene_ensembl_file = file(params.gene_ensembl_file)
            assert gene_ensembl_file.exists() : "Error: gene_ensembl_file not found: ${params.gene_ensembl_file}"
            
            log.info "📊 Generating CAAS characterization reports..."
            
            // Pass the filter_ch output directory path instead of individual files
            def filter_output_dir = params.caas_postproc_mode == 'exploratory' ? 
                "${params.outdir}/postproc/filter_${params.caas_postproc_mode}" :
                "${params.outdir}/postproc/filter_selected"
            
            characterization_results = CAAS_POSTPROC_REPORT(
                discovery_file_ch,
                filter_summary.summary,
                filter_output_dir,
                gene_ensembl_file
            )
            
            log.info "Post-processing reports generated in: ${params.outdir}/postproc/reports"
        }
    
    emit:
        filter_summary = filter_summary.summary
        filter_dir = "${params.outdir}/postproc"
        filtered_discovery = gene_filter_results ? gene_filter_results.filtered_discovery : Channel.empty()
        cleaned_background = cleaned_backgrounds
}
