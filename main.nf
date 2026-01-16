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
    } else if (params.ct_tool){
        CT ()
    } else if (params.rer_tool) {
        RER_MAIN()
    } /* else if (params.ora) {
        ORA()
    } */
}


/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  THE END: End of the main.nf file.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
