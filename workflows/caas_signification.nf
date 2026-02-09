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
    
    main:
        // Input validation
        assert params.bootstrap_input : "CAAS Signification requires --bootstrap_input. Provide path to bootstrap file or directory"
        
        def bootstrap_path = file(params.bootstrap_input)
        assert bootstrap_path.exists() : "Error: bootstrap_input not found: ${params.bootstrap_input}"
        
        // Handle bootstrap files: single file or directory
        def bootstrap_files
        def bootstrap_count
        
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
        
        // Detect CAAP mode from filtered discovery file
        def is_caap_mode = params.caap_mode ?: false
        
        // Log configuration
        log.info """
        ╔════════════════════════════════════════════════════════════════╗
        ║            CAAS SIGNIFICATION WORKFLOW                         ║
        ╚════════════════════════════════════════════════════════════════╝
        
        CAAP Mode:            ${is_caap_mode ? 'Yes (group-aware analysis)' : 'No (standard analysis)'}
        Bootstrap Files:      ${bootstrap_count} files from ${params.bootstrap_input}
        Output Directory:     ${params.signification_outdir}
        Significance Alpha:   0.05 (hypergeometric & permutation)
        
        Analysis Steps:
          1. Load filtered CAAS discovery data
          2. Integrate bootstrap permutation p-values
          3. Apply significance filters (hypergeometric, permutation, both)
          4. Analyze pattern distributions
          5. Export gene lists for enrichment
          6. Generate meta-CAAS metadata
        
        """.stripIndent()
        
        // Run signification analysis
        CAAS_SIGNIFICATION_REPORT(
            filtered_discovery_channel,
            filtered_background_channel,
            bootstrap_files
        )
        
        // Log completion
        log.info "✓ CAAS Signification analysis complete"
        log.info "  Reports available in: ${params.signification_outdir}"
        
    emit:
        report = CAAS_SIGNIFICATION_REPORT.out.report
        gene_lists = CAAS_SIGNIFICATION_REPORT.out.gene_lists
        commonalities = CAAS_SIGNIFICATION_REPORT.out.commonalities
}
