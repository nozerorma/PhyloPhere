#!/usr/bin/env nextflow

/*
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#                                                                                    
#                                      
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: main.nf
#
*/

/*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Unlock the secrets of evolutionary relationships with Phylophere! 🌳🔍 This Nextflow pipeline
* packs a powerful punch, offering a comprehensive suite of phylogenetic comparative tools
* and analyses. Dive into the world of evolutionary biology like never before and elevate
* your research to new heights! 🚀🧬 #Phylophere #EvolutionaryInsights #NextflowPipeline
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

version = "0.0.1"

// Display input parameters
log.info """

PHYLOPHERE - NF PIPELINE  ~  version ${version}
=============================================

PHYLOPHERE: A Nextflow pipeline including a complete set
of phylogenetic comparative tools and analyses for Phenome-Genome studies

Author:         Miguel Ramon (miguel.ramon@upf.edu)


"""

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  NAMED WORKFLOW FOR PIPELINE: This section includes the main CT, ORA and RERConverge workflows.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include {HELP} from './workflows/help.nf'
include {CT} from './workflows/ct.nf'
include {RER_MAIN} from './workflows/rerconverge.nf'
include {REPORTING} from './workflows/reporting.nf'
include {CONTRAST_SELECTION} from './workflows/contrast_selection.nf'
include {CAAS_POSTPROC} from './workflows/caas_postproc.nf'
include {CAAS_SIGNIFICATION} from './workflows/caas_signification.nf'
//include {ORA} from './workflows/ora.nf'

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RUN PHYLOPHERE ANALYSIS: This section initiates the main Phylophere workflow.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow {

    // Check if --help is provided
    if (params.help) {
        HELP ()
    } else {
        // Run any combination of tools requested
        def ran_any = false

        if (params.reporting && !params.contrast_selection) {
            REPORTING()
            ran_any = true
        }
        def ct_results
        if (params.ct_tool) {
            if (params.contrast_selection) {
                def contrast_out = CONTRAST_SELECTION()
                ct_results = CT(contrast_out.trait_file_out, contrast_out.bootstrap_trait_file_out, contrast_out.tree_file_out)
            } else {
                def trait_file_in = null
                def bootstrap_trait_file_in = null
                def tree_file_in = null
                ct_results = CT (trait_file_in, bootstrap_trait_file_in, tree_file_in)
            }
            ran_any = true
        }
        if (params.contrast_selection && !params.ct_tool) {
            CONTRAST_SELECTION()
            ran_any = true
        }
        if (params.rer_tool) {
            RER_MAIN()
            ran_any = true
        }
        if (params.caas_postproc) {
            // Use CT outputs if available, otherwise use empty channels (will fall back to parameters)
            def discovery_ch = ct_results ? ct_results.discovery_file : Channel.empty()
            def background_ch = ct_results ? ct_results.background_file : Channel.empty()
            def background_genes_ch = ct_results ? ct_results.background_genes : Channel.empty()
            def bootstrap_ch = ct_results ? ct_results.bootstrap_file : Channel.empty()
            
            def postproc_results = CAAS_POSTPROC(discovery_ch, background_ch, background_genes_ch)
            ran_any = true
            
            // Optionally chain CAAS_SIGNIFICATION after CAAS_POSTPROC
            if (params.caas_signification) {
                // For signification, we need bootstrap data from CT if available
                // Otherwise it will use params.bootstrap_input
                CAAS_SIGNIFICATION(
                    postproc_results.filtered_discovery,
                    postproc_results.cleaned_background,
                    bootstrap_ch
                )
            }
        }
        if (params.caas_signification && !params.caas_postproc) {
            log.error "ERROR: --caas_signification requires --caas_postproc to run first"
            log.error "Please enable both workflows or run CAAS_POSTPROC separately to generate filtered inputs"
            System.exit(1)
        }
        /* if (params.ora) {
            ORA()
            ran_any = true
        } */

        if (!ran_any) {
            log.info "No tool selected. Use --reporting, --contrast_selection, --ct_tool, --rer_tool, --caas_postproc, or --caas_signification."
        }
    }
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  THE END: End of the main.nf file.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
