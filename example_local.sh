#!/bin/bash
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
# File: discovery-local.sh
#

set -Eeuo pipefail

timestamp=$(date +%Y%m%d%H%M%S)

isToy=TRUE
traits_to_run=("neoplasia_prevalence" "malignant_prevalence")

# Alignment dir in the specified format. Toy is just subset.
if [ $isToy = TRUE ]; then
    tag="_toy"
    ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Ali_toy"
    CYCLES="100"
else
    tag=""
    ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Primate_alignments"
    CYCLES="1000000"
fi


# ## CAASTOOLS DISCOVERY
# TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Traitfiles/$TRAIT/traitfile.tab" # Directory where your trait files are located

# ## CAASTOOLS RESAMPLE
# TREE_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree" # Path to the tree file
# TRAIT_VALUES="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.5.Bootstrap_traitfiles/$TRAIT/boot_traitfile.tab" # Path to the trait values file
# CHUNK_SIZE="500" # Cycles per file (creates directory with multiple resample_*.tab files)

## CAASTOOLS BOOTSTRAP
#RESAMPLED_DIR="/path/goes/here/user/Resampled_files"                # If resample already done, path to the resampled directory
#DISCOVERY_DIR="/path/to/discovery/output"                           # Optional: for position filtering optimization (100-1000x speedup)
#PROGRESS_LOG="bootstrap_progress.log"                                # Optional: progress tracking with timestamps and ETA


## If running locally (not using Singularity)
CT_PATH="./caastools/ct"                                      # Path to the CAASTOOLS binary

# Pipeline execution

### If using Singularity, add the -with-singularity flag and the -profile singularity flag
### If using Nextflow Tower, add the -with-tower flag

# Apptainer run
## nextflow run main.nf -with-tower -with-apptainer -profile local,apptainer

# Local run
## nextflow run main.nf -with-tower -profile local

# Slurm localized run (using haswell, with apptainer)
## nextflow run main.nf -with-tower -with-apptainer -profile local,apptainer

# Slurm delocalized run (using haswell, using slurm, with apptainer). Careful with slurm queueing, not very efficient due to overhead.
# Queue should be set in the sbatch directly.
## nextflow run main.nf -with-tower -with-apptainer -profile slurm,apptainer

# Singularity is legacy. Apptainer should be used unless necessary.
# For some reason, Nextflow does not mount HOME directly no more. We need to flag it.
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true


for trait in "${traits_to_run[@]}"; do
    # CAASTOOLS DISCOVERY
    echo "Running PhyloPhere for trait: $trait"
    RESULTS_DIR="Out/CAAS_2.0/${trait}${tag}"   # Directory where results will be stored
    mkdir -p $RESULTS_DIR


    WORK_DIR="Out/CAAS_2.0/${trait}${tag}/Work"
    mkdir -p $WORK_DIR
    
    out_report="$RESULTS_DIR/report/$timestamp"
    mkdir -p $out_report

    ## CAASTOOLS DISCOVERY
    TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Traitfiles/$trait/traitfile.tab" # Directory where your trait files are located

    ## CAASTOOLS RESAMPLE
    TREE_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree" # Path to the tree file
    TRAIT_VALUES="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.5.Bootstrap_traitfiles/$trait/boot_traitfile.tab" # Path to the trait values file
    
    # nextflow run main.nf -with-tower -profile local \
    #     -w $WORK_DIR \
    #     --reporting \
    #     --outdir $out_report

    # sleep 2s

    out_analysis="$RESULTS_DIR/analysis_prune/$timestamp"
    mkdir -p $out_analysis

    nextflow run main.nf -with-tower -profile local \
        -w $WORK_DIR \
        --prune_data \
        --contrast_selection \
        --ct_tool "discovery,resample,bootstrap" \
        --alignment $ALI_DIR \
        --caas_config $TRAIT_FILE \
        --tree $TREE_FILE \
        --traitvalues $TRAIT_VALUES \
        --cycles $CYCLES \
        --outdir $out_analysis 

    out_analysis="$RESULTS_DIR/analysis_raw/$timestamp"
    mkdir -p $out_analysis

    nextflow run main.nf -with-tower -profile local \
        -w $WORK_DIR \
        --contrast_selection \
        --ct_tool "discovery,resample,bootstrap" \
        --alignment $ALI_DIR \
        --caas_config $TRAIT_FILE \
        --tree $TREE_FILE \
        --traitvalues $TRAIT_VALUES \
        --cycles $CYCLES \
        --outdir $out_analysis 

    out_analysis="$RESULTS_DIR/analysis_file/$timestamp"
    mkdir -p $out_analysis

    nextflow run main.nf -with-tower -profile local \
        -w $WORK_DIR \
        --ct_tool "discovery,resample,bootstrap" \
        --alignment $ALI_DIR \
        --caas_config $TRAIT_FILE \
        --tree $TREE_FILE \
        --traitvalues $TRAIT_VALUES \
        --cycles $CYCLES \
        --outdir $out_analysis 

    echo "Phylophere finished the analysis for trait: $trait."
    echo "Results stored in: $RESULTS_DIR"
done

#debug flags
# --export_groups  # Export the groups used in the bootstrap resample for debugging
# --export_perm_discovery  # Export the discovery results for each permutation in the resample for debugging
# --include_b0 # Export b_0 (hypothesis test) in the resample for debugging
# --resample_out # Path where resampled files are stored for independent bootstrap runs
# --discovery_out # Path where discovery output files are stored for independent bootstrap runs
#    --resample_out "/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/CAAS_2.0/1.Discovery/neoplasia_prevalence_toy/20260122183425/resample/science.abn7829_data_s4.nex.resampled.output"

# nextflow clean -f # This flag should be disabled if debugging
echo "---------------------------------------------"
echo "Pipeline finished successfully."
echo "---------------------------------------------"
echo "Enjoy!"
echo "End of analysis at $(timestamp)"
