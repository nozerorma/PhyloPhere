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
        background_files_channel  // Background files from CT module (required if not using --background_input)
    
    main:
        // Input validation
        assert params.discovery_input : "CAAS Post-Processing requires --discovery_input. Provide path to discovery.tab file"
        
        def discovery_file = file(params.discovery_input)
        
        assert discovery_file.exists() : "Error: discovery_input file not found: ${params.discovery_input}"
        assert discovery_file.isFile() : "Error: discovery_input must be a file"
        
        // Detect CAAP mode by checking for CAAP_Group column in discovery file
        def is_caap_mode = false
        def caap_groups = []
        
        // Read first few lines to check for CAAP_Group column
        discovery_file.withReader { reader ->
            def header = reader.readLine()
            if (header && header.contains("CAAP_Group")) {
                is_caap_mode = true
                log.info "📂 CAAP mode detected: discovery.tab contains CAAP_Group column"
            } else {
                log.info "📄 CAAS mode detected: standard discovery.tab format"
            }
        }
        
        // Handle background files: use parameter if provided, otherwise use channel input
        def background_files
        if (params.background_input) {
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
            // Use channel input from CT module (must not be empty)
            background_files = background_files_channel
            log.info "📥 Using background files from CT_DISCOVERY module"
        }
        
        // Log configuration
        log.info """
        ╔════════════════════════════════════════════════════════════════╗
        ║            CAAS POST-PROCESSING WORKFLOW                       ║
        ╚════════════════════════════════════════════════════════════════╝
        
        Mode:                 ${params.caas_postproc_mode}
        CAAP Mode:            ${is_caap_mode ? 'Yes (multiple groups)' : 'No (single file)'}
        Discovery Input:      ${params.discovery_input}
        Gene Ensembl File:    ${params.gene_ensembl_file ?: 'N/A (reports disabled)'}
        Background Input:     ${params.background_input ?: 'From CT_DISCOVERY channel (required)'}
        Gene Filter Mode:     ${params.gene_filter_mode}
        Output Directory:     ${params.postproc_outdir}
        Generate Reports:     ${params.generate_reports}
        Generate Manhattan:   ${params.generate_manhattan}
        
        """.stripIndent()
        
        // Determine processing mode and create parameter combinations channel
        if (params.caas_postproc_mode == 'exploratory') {
            // Parameter sweep: generate cartesian product of all parameter combinations
            def minlen_list = params.minlen_values.split(',').collect { it.trim().toInteger() }
            def maxcaas_list = params.maxcaas_values.split(',').collect { it.trim().toDouble() }
            
            // Create channel with all parameter combinations
            param_combinations = Channel
                .from(minlen_list)
                .combine(Channel.from(maxcaas_list))
                .map { minlen, maxcaas -> 
                    tuple('exploratory', minlen, maxcaas, discovery_file)
                }
            
            log.info "🔍 Exploratory mode: testing ${minlen_list.size()} × ${maxcaas_list.size()} = ${minlen_list.size() * maxcaas_list.size()} parameter combinations"
            if (is_caap_mode) {
                log.info "   (Cluster filtering will respect CAAP_Group boundaries)"
            }
            
        } else if (params.caas_postproc_mode == 'filter') {
            // Single filter run: use provided minlen and maxcaas
            param_combinations = Channel.of(
                tuple('filter', params.filter_minlen, params.filter_maxcaas, discovery_file)
            )
            
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
            filter_results.filtered.map { minlen, maxcaas, file -> file }.collect()
        )
        
        // Run gene-level filtering if enabled
        def gene_filter_results
        def cleaned_backgrounds
        
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
            def cluster_file = filter_results.filtered.map { minlen, maxcaas, file -> file }.first()
            
            gene_filter_results = CAAS_FILTER_GENES(
                discovery_file,
                gene_ensembl_file,
                cluster_file
            )
            
            // Run background cleanup (always when gene filtering is enabled)
            log.info "🧹 Cleaning background gene lists..."
            
            cleaned_backgrounds = CAAS_BACKGROUND_CLEANUP(
                background_files,
                gene_filter_results.removed_genes
            )
            
            log.info "✓ Gene filtering complete"
            log.info "✓ Cleaned background files: ${params.postproc_outdir}/cleaned_backgrounds"
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
            
            // Single discovery file with group-aware reporting
            characterization_results = CAAS_POSTPROC_REPORT(
                discovery_file,
                filter_summary.summary,
                gene_ensembl_file
            )
            
            log.info "✓ Reports generated in: ${params.postproc_outdir}/reports"
        }
    
    emit:
        filter_summary = filter_summary.summary
        filter_dir = params.postproc_outdir
        filtered_discovery = gene_filter_results ? gene_filter_results.filtered_discovery : Channel.empty()
        cleaned_background = cleaned_backgrounds
}
