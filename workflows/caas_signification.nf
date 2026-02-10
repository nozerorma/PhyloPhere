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
# File: caas_signification.nf
#
*/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  CAAS Signification Workflow: Performs significance testing using hypergeometric and
 *  permutation tests with bootstrap data. Exports gene lists for enrichment analysis
 *  and meta-CAAS data for cross-study comparison.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Import local processes from subworkflows
include { CAAS_SIGNIFICATION_REPORT } from "${baseDir}/subworkflows/CAAS_SIGNIFICATION/ctpp_signification"

workflow CAAS_SIGNIFICATION {
    take:
        filtered_discovery_channel  // Filtered discovery file from CAAS_POSTPROC
        filtered_background_channel // Cleaned background file from CAAS_POSTPROC
        bootstrap_input_channel     // Bootstrap file from CT module (optional, can use --bootstrap_input instead)
    
    main:
        // Determine bootstrap file source: channel input or parameter
        def bootstrap_files
        def bootstrap_count
        
        if (params.ct_tool && bootstrap_input_channel) {
            log.info "📥 Using bootstrap file from CT module"
            bootstrap_files = bootstrap_input_channel.map { file -> file }.collect()
            bootstrap_count = 1
        } else {
            assert params.bootstrap_input : "CAAS Signification requires --bootstrap_input when not using CT pipeline. Provide path to bootstrap file or directory"
            
            def bootstrap_path = file(params.bootstrap_input)
            assert bootstrap_path.exists() : "Error: bootstrap_input not found: ${params.bootstrap_input}"
        
            // Handle bootstrap files: single file or directory
            if (bootstrap_path.isDirectory()) {
                // Directory: collect all .boot and .tab files
                def boot_files = Channel.fromPath("${params.bootstrap_input}/*.{boot,tab}")
                bootstrap_files = boot_files.collect()
                bootstrap_count = file(params.bootstrap_input).list().findAll { 
                    it.endsWith('.boot') || it.endsWith('.tab') 
                }.size()
                assert bootstrap_count > 0 : "Error: No .boot or .tab files found in directory ${params.bootstrap_input}"
                log.info "📂 Loading ${bootstrap_count} bootstrap files from directory: ${params.bootstrap_input}"
            } else {
                // Single file (accept .boot or .tab extension)
                assert bootstrap_path.name.endsWith('.boot') || bootstrap_path.name.endsWith('.tab') : \
                    "Error: bootstrap_input file must have .boot or .tab extension"
                bootstrap_files = Channel.fromPath(params.bootstrap_input).collect()
                bootstrap_count = 1
                log.info "📄 Loading single bootstrap file: ${params.bootstrap_input}"
            }
        }
        
        // Detect CAAP mode from filtered discovery file
        def is_caap_mode = params.caap_mode ?: false
        
        // Run signification analysis
        CAAS_SIGNIFICATION_REPORT(
            filtered_discovery_channel,
            filtered_background_channel,
            bootstrap_files
        )
        
        // Log completion
        log.info "Signification reports available in: ${params.outdir}/signification"
        
    emit:
        report = CAAS_SIGNIFICATION_REPORT.out.report
        gene_lists = CAAS_SIGNIFICATION_REPORT.out.gene_lists
        meta_caas = CAAS_SIGNIFICATION_REPORT.out.meta_caas
}
