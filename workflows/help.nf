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
# File: help.nf
#

#Detects Candidate Amino Acid Substitutions (CAAS) from a single Multiple Sequence Alignment (MSA).

Parameters:
--alignment               Path to the alignment file. Default: "$baseDir/examples/MSA"
--disc_description        Description for the discovery. Default: null
--traitfile               Path to the trait file. Default: "$baseDir/examples/config.tab"
... [rest of the shared parameters between discovery and bootstrap]
*/

// General help message
def general_help = '''
CAASTOOLS Nextflow pipeline
=============================================
Convergent Amino Acid Substitution detection and analysis TOOLbox
in a Nextflow fashion.

General usage:          > nextflow run main.nf [options]
Help for single tool:   > nextflow run main.nf --ct_tool <tool> --help

Tools           Description
--------        -------------------------------------------------
discovery       Detects Candidate Amino Acid Substitutions (CAAS) from
                a single Multiple Sequence Alignment (MSA).

resample        Resamples virtual phenotypes for CAAS bootstrap analysis.

bootstrap       Runs CAAS bootstrap analysis on a on a single MSA.

pgls            Runs site-level phylogenetic GLS on CT post-processed CAAS
                sites and renders an HTML report.

scoring         Computes composite CAAS scores at position-level and
                gene-level, integrating CT, PGLS, FADE, RER, and
                accumulation signals.

'''

// Define tool-specific help messages
def discovery_help = '''
Discovery Tool Help
=============================================
Detects Candidate Amino Acid Substitutions (CAAS) from a single Multiple Sequence Alignment (MSA).

Usage:
--alignment             <"input_dir">                           null
--traitfile             <"traitfile">                           null
--ali_format            <"ali_format">                          null
--patterns              <"1,2,3,4">                             "1,2,3"
--maxbggaps             <"NO|INTEGER">                          "NO"
--maxfggaps             <"NO|INTEGER">                          "NO"
--maxgaps               <"NO|INTEGER">                          "NO"
--maxgapsperposition    <"INTEGER">                             "0.5"
--maxbgmiss             <"NO|INTEGER">                          "NO"
--maxfgmiss             <"NO|INTEGER">                          "NO"
--maxmiss               <"NO|INTEGER">                          "NO" 
'''

def resample_help = '''
Resample Tool Help
=============================================
Resamples virtual phenotypes for CAAS bootstrap analysis.

Usage:
--tree                  <"nwtree_file">                         ${params.tree}
--strategy              <"FGBG|BM">                             ${params.strategy}
--fgsize                <"INTEGER">                             ${params.fgsize}
--bgsize                <"INTEGER">                             null
--traitvalues           <"traitvalues_file">                    null
--cycles                <"INTEGER">                             "1000"
--chunk_size            <"INTEGER">                             "500"
--include_b0            <"true|false">                          true

Output: Directory containing resample_*.tab files (one per chunk_size cycles)

Strategy requirements:
FGBG                    --fgsize --bgsize
BM                      --bytemp --traitvalues
'''

def bootstrap_help = '''
Bootstrap Tool Help
=============================================
Runs CAAS bootstrap analysis on a single MSA.
Usage:
--resample_out          <"resampleDir|resampleFile">            null
--progress_log          <"progress_log_file">                   "none"
--discovery_out         <"discovery_output_dir">                "none"
--export_groups         <"output_file">                         "none"
--export_perm_discovery <"output_file">                         "none"

NOTE: resample_out can be either a directory (recommended, contains resample_*.tab files) or a legacy single file
NOTE: discovery_out enables position filtering optimization (typically 100-1000× speedup)
NOTE: progress_log creates timestamped progress tracking with ETA

# Common parameters with alignment tool
--alignment             <"input_dir">                           null
--traitfile             <"traitfile">                           null
--ali_format            <"ali_format">                          null
--patterns              <"1,2,3,4">                             "1,2,3"
--maxbggaps             <"NO|INTEGER">                          "NO"
--maxfggaps             <"NO|INTEGER">                          "NO"
--maxgaps               <"NO|INTEGER">                          "NO"
--maxgapsperposition    <"INTEGER">                             "0.5"
--maxbgmiss             <"NO|INTEGER">                          "NO"
--maxfgmiss             <"NO|INTEGER">                          "NO"
--maxmiss               <"NO|INTEGER">                          "NO" 
'''

def pgls_help = """
Site-PGLS Help
=============================================
Runs site-level phylogenetic GLS association tests on CT post-processed CAAS sites.

Usage:
--pgls                  <true|false>                            true
--pgls_caas_input       <"filtered_discovery.tsv">              null
--my_traits             <"traits.csv">                          ${params.my_traits}
--traitname             <"trait_column">                        ${params.traitname}
--tree                  <"tree_file">                           ${params.tree}
--alignment             <"alignment_dir">                       ${params.alignment}
--ali_format            <"alignment_format">                    ${params.ali_format}
--tax_id                <"tax_map.tsv">                         ${params.tax_id}
--sp_colname            <"species_column">                      ${params.sp_colname}
--discrete_method       <"quartile|quintile|decile|median_sd|parameterized"> ${params.discrete_method}
--top_quantile          <"FLOAT">                               ${params.top_quantile}
--bottom_quantile       <"FLOAT">                               ${params.bottom_quantile}
--use_median_sep        <"true|false">                          ${params.use_median_sep}
--pgls_threads          <"INTEGER">                             ${params.pgls_threads}
--pgls_site_fdr_alpha   <"FLOAT">                               ${params.pgls_site_fdr_alpha}
--pgls_min_per_class    <"INTEGER">                             ${params.pgls_min_per_class}
--pgls_enrichment       <"true|false">                          ${params.pgls_enrichment}
--pgls_phylo_model      <"bm|lambda|ou|eb|delta|kappa">         ${params.pgls_phylo_model}
"""
def scoring_help = """
CAAS Scoring Help
=============================================
Computes composite CAAS scores at position-level and gene-level.

Usage:
--scoring                   <true|false>            false
--scoring_ora               <true|false>            true
--scoring_stress            <true|false>            false
--scoring_stress_top_n      <INTEGER>               25
--scoring_stress_rank_metric <"spearman">           spearman
--scoring_position_top_pct  <FLOAT>                 0.10
--scoring_position_top5_pct <FLOAT>                 0.05
--scoring_position_top1_pct <FLOAT>                 0.01
--scoring_gene_top_pct      <FLOAT>                 0.10
--scoring_gene_top5_pct     <FLOAT>                 0.05
--scoring_gene_top1_pct     <FLOAT>                 0.01

Standalone mode (provide inputs directly):
--scoring_postproc_input    <"filtered_discovery.tsv">   ""
--scoring_pgls_input        <"site_pgls.tsv">            ""
--scoring_pgls_excess_input <"pgls_excess_gene_trait.tsv"> ""
--scoring_fade_summary_top  <"fade_summary_top.tsv">     ""
--scoring_fade_summary_bottom <"fade_summary_bottom.tsv"> ""
--scoring_rer_input         <"rerconverge_summary.tsv">  ""
--scoring_accum_dir                 <"accumulation_dir/">              ""
--scoring_molerate_summary_top      <"molerate_summary_top.tsv">       ""
--scoring_molerate_summary_bottom   <"molerate_summary_bottom.tsv">    ""
--scoring_background_input          <"background.txt">                 ""

Position-level components: biochem, ASR, convergence, parallel, [PGLS], [FADE]
Gene-level scores: gene_caas, [gene_rand], [gene_rer], [gene_fade], [gene_molerate], gene_composite
"""

workflow HELP {
    // Check if --help is provided
    if (params.help) {
        if (params.scoring) {
            log.info scoring_help
            exit 1
        }
        if (params.pgls) {
            log.info pgls_help
            exit 1
        }
        // Check if a specific tool is mentioned with --ct_tool
        if (params.ct_tool) {
            switch (params.ct_tool) {
                case 'discovery':
                    log.info discovery_help
                    break
                case 'resample':
                    log.info resample_help
                    break
                case 'bootstrap':
                    log.info bootstrap_help
                    break
                default:
                    log.info general_help
            }
        } else {
            // If no specific tool is mentioned, display the general help message
            log.info general_help
        }
        exit 1
    }
}
