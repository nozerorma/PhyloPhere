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

version = "1.0.0"

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
include {CT_SIGNIFICATION} from './workflows/ct_signification.nf'
include {CT_POSTPROC} from './workflows/ct_postproc.nf'
include {CT_DISAMBIGUATION} from './workflows/ct_disambiguation.nf'
include {ORA} from './workflows/ora.nf'
include {ORA_ACCUMULATION} from './workflows/ora_accumulation.nf'
include {CT_ACCUMULATION} from './workflows/ct_accumulation.nf'
include {FADE}           from './workflows/fade.nf'
include {MOLERATE}       from './workflows/molerate.nf'

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

                // Hard stop: if CHECK_MIN_CONTRASTS emits low_contrasts.skip,
                // terminate the current trait run gracefully (exit 0).
                contrast_out.low_contrasts_skip.view { skip_file ->
                    exit 0, "Minimum contrast threshold not met for trait '${params.traitname ?: 'unknown'}' (flag: ${skip_file}). Stopping pipeline gracefully."
                }

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
        def signification_results = null
        def disambiguation_results = null
        def postproc_results = null

        // Channels populated by CT_ACCUMULATION/CT_POSTPROC when they run.
        // Used by FADE/MOLERATE/RER gene-set piping without manual path specification.
        def sel_acc_top_ch    = Channel.empty()
        def sel_acc_bottom_ch = Channel.empty()
        def sel_pp_top_ch     = Channel.empty()
        def sel_pp_bottom_ch  = Channel.empty()

        if (params.ct_signification) {
            // Use CT outputs if available, otherwise workflow falls back to params inputs.
            def discovery_ch       = ct_results ? ct_results.discovery_file  : null
            def background_genes_ch = ct_results ? ct_results.background_genes : null
            def bootstrap_ch       = ct_results ? ct_results.bootstrap_file   : null

            signification_results = CT_SIGNIFICATION(discovery_ch, background_genes_ch, bootstrap_ch)
            ran_any = true
        }

        if (params.ct_disambiguation) {
            if (!params.ct_signification && !params.ct_disambig_caas_metadata) {
                error "CT disambiguation requires signification outputs (--ct_signification) or a standalone metadata file (--ct_disambig_caas_metadata)."
            }

            // Pass null so that the if(channel) guard inside CT_DISAMBIGUATION detects absence
            // and falls back to --ct_disambig_caas_metadata.
            def meta_for_disambiguation = signification_results ? signification_results.signification_global_meta : null
            def trait_for_disambiguation = ct_results ? ct_results.trait_file : Channel.empty()
            def tree_for_disambiguation = ct_results ? ct_results.tree_file : Channel.empty()

            disambiguation_results = CT_DISAMBIGUATION(meta_for_disambiguation, trait_for_disambiguation, tree_for_disambiguation)
            ran_any = true
        }

        if (params.ora && !params.ct_postproc) {
            // Standalone ORA: no upstream postproc channels — ORA will fall back to
            // --ora_gene_lists_input + --ora_background_input via its .ifEmpty {} guards.
            if (!params.ora_gene_lists_input || !params.ora_background_input) {
                error "Standalone ORA requires --ora_gene_lists_input and --ora_background_input (or enable --ct_postproc for integrated mode)."
            }
            // Pass null for gene_lists so the if(channel) guard in ora.nf detects absence
            // and loads directly from --ora_gene_lists_input (Channel.empty() would silently
            // stall .collect() without ever emitting).
            ORA(Channel.empty(), null)
            ran_any = true
        }

        if (params.ct_postproc) {
            if (!params.ct_disambiguation && !params.discovery_input) {
                error "CT post-processing now runs downstream of disambiguation. Enable --ct_disambiguation or provide --discovery_input (caas_convergence_master.csv)."
            }

            // Post-processing is downstream from disambiguation; consume disambiguation master CSV when available
            // Pass null (not Channel.empty()) when there is no upstream result so that the
            // if(channel) guard inside CT_POSTPROC correctly detects absence and falls back
            // to --discovery_input / --background_input params (same pattern as CT_SIGNIFICATION).
            def discovery_ch = disambiguation_results ? disambiguation_results.master_csv : null
            def background_ch = ct_results ? ct_results.background_file_raw : Channel.empty()
            def background_genes_ch = ct_results ? ct_results.background_genes : null
            def bootstrap_ch = Channel.empty() // retained for CT_POSTPROC signature compatibility

            postproc_results = CT_POSTPROC(discovery_ch, background_ch, background_genes_ch, bootstrap_ch)
            ran_any = true

            if (params.ora) {
                // ORA consumes CT_POSTPROC channels in integrated runs
                def ora_background_ch = postproc_results.cleaned_background
                def ora_gene_lists_ch = postproc_results.ora_gene_lists_files

                ORA(ora_background_ch, ora_gene_lists_ch)
                ran_any = true
            }
        }

        if (params.ct_accumulation) {
            if (!params.ct_postproc && !params.accumulation_background_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone background file (--accumulation_background_input)."
            }
            if (!params.ct_postproc && !params.accumulation_caas_input) {
                error "CT_ACCUMULATION requires CT post-processing output (--ct_postproc) or a standalone CAAS file (--accumulation_caas_input)."
            }

            // Use filtered_discovery.tsv from postproc (gene_filtering stage)
            def acc_caas_ch       = postproc_results      ? postproc_results.filtered_discovery             : Channel.empty()
            def acc_background_ch = postproc_results      ? postproc_results.cleaned_background             : Channel.empty()
            def acc_trait_file_ch = ct_results            ? ct_results.trait_file                          : Channel.empty()

            def accum_results = CT_ACCUMULATION(acc_caas_ch, acc_background_ch, acc_trait_file_ch)
            ran_any = true

            // Split accumulation CSVs by direction for FADE/MOLERATE gene-set piping
            sel_acc_top_ch    = accum_results.results
                .filter { f -> f.name.contains('_top_') }
            sel_acc_bottom_ch = accum_results.results
                .filter { f -> f.name.contains('_bottom_') }

            if (params.ora) {
                // Pipe accumulation gene lists into the ORA / STRING enrichment pipeline.
                // Background falls back to --ora_background_input when no postproc channel is present.
                def ora_acc_background_ch  = postproc_results ? postproc_results.cleaned_background : Channel.empty()
                ORA_ACCUMULATION(ora_acc_background_ch, accum_results.gene_lists)
            }
        }

        // Populate postproc gene-list channels for FADE/MOLERATE gene-set piping
        if (postproc_results) {
            sel_pp_top_ch    = postproc_results.ora_gene_lists_files
                .filter { f -> f.name.contains('_top_') && f.name.endsWith('significant.txt') }
            sel_pp_bottom_ch = postproc_results.ora_gene_lists_files
                .filter { f -> f.name.contains('_bottom_') && f.name.endsWith('significant.txt') }
        }

        if (params.fade) {
            def fade_traitfile_ch = ct_results ? ct_results.trait_file : Channel.empty()
            def fade_tree_ch      = ct_results ? ct_results.tree_file  : Channel.empty()
            FADE(fade_traitfile_ch, fade_tree_ch,
                 sel_acc_top_ch, sel_acc_bottom_ch,
                 sel_pp_top_ch,  sel_pp_bottom_ch)
            ran_any = true
        }

        if (params.molerate) {
            def mr_traitfile_ch = ct_results ? ct_results.trait_file : Channel.empty()
            def mr_tree_ch      = ct_results ? ct_results.tree_file  : Channel.empty()
            MOLERATE(mr_traitfile_ch, mr_tree_ch,
                     sel_acc_top_ch, sel_acc_bottom_ch,
                     sel_pp_top_ch,  sel_pp_bottom_ch)
            ran_any = true
        }

        if (params.rer_tool) {
            // RER_MAIN is called last so that sel_acc_*/sel_pp_* channels are
            // fully populated by CT_ACCUMULATION / CT_POSTPROC when running together.
            def rer_traitfile_ch = ct_results ? ct_results.trait_file : Channel.empty()
            RER_MAIN(
                rer_traitfile_ch,
                sel_acc_top_ch, sel_acc_bottom_ch,
                sel_pp_top_ch,  sel_pp_bottom_ch
            )
            ran_any = true
        }

        if (!ran_any) {
            log.info "No tool selected. Use --reporting, --contrast_selection, --ct_tool, --rer_tool, --ct_signification, --ct_disambiguation, --ct_postproc, --ora, --ct_accumulation, --fade, --molerate, or --rer_tool."
        }
    }
}

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  THE END: End of the main.nf file.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
