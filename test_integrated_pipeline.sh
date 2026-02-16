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
# PHYLOPHERE: Integrated Pipeline Test Script
#
# This script runs the complete pipeline including:
# 1. CAAS Discovery (CT tool)
# 2. Resample generation
# 3. Bootstrap analysis
# 4. CAAS Post-Processing (filtering and characterization)
#
# The outputs from steps 1-3 are automatically piped into step 4,
# demonstrating the full end-to-end workflow.
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: test_integrated_pipeline.sh
#

set -Eeuo pipefail

# Activate micromamba environment
echo "Activating micromamba environment: caas_global_cancer"
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate caas_global_cancer
else
    echo "Warning: micromamba not found in PATH, attempting manual conda activation"
    if command -v conda &> /dev/null; then
        conda activate caas_global_cancer || echo "Warning: Could not activate conda environment"
    else
        echo "Warning: Neither micromamba nor conda found, proceeding without environment activation"
    fi
fi

timestamp=$(date +%Y%m%d_%H%M%S)

# Configuration
# Set to true to run a quick test with a small subset of data (for debugging purposes)
IS_TOY=true

# Base data directory
DATADIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data"

# Trait and branch definitions
TRAIT="neoplasia_prevalence"
SECONDARY_TRAIT="malignant_prevalence"
N_TRAIT="adult_necropsy_count"
C_TRAIT="neoplasia_necropsy"
BRANCH_TRAIT="LQ"

# Pruning lists
PRUNE_LIST="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP/neoplasia_exclude.txt"
PRUNE_SECONDARY_LIST="${DATADIR}/1.Cancer_data/Neoplasia_species360/ZAK-CLEANUP/malignant_exclude.txt"

# Trait file
TRAIT_FILE="${DATADIR}/1.Cancer_data/Neoplasia_species360/cancer_traits_processed-LQ.csv"


# Input configuration
if [ "$IS_TOY" = true ]; then
    TAG="_toy"
    ALI_DIR="${DATADIR}/2.Alignments/Ali_toy"
    CYCLES="100"
else
    TAG=""
    ALI_DIR="${DATADIR}/2.Alignments/Primate_alignments"
    CYCLES="1000000"
fi


# Validate required files
for file in "$TRAIT_FILE"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file not found: $file"
        exit 1
    fi
done

if [ ! -d "$ALI_DIR" ]; then
    echo "Error: Alignment directory not found: $ALI_DIR"
    exit 1
fi

# Output configuration
RESULTS_BASE="/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/CAAS_2.0/${TRAIT}${TAG}/${timestamp}"
mkdir -p "$RESULTS_BASE"

WORK_DIR="${RESULTS_BASE}/work"
mkdir -p "$WORK_DIR"

# For Nextflow Tower/Apptainer
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true

echo "=========================================="
echo "PHYLOPHERE INTEGRATED PIPELINE TEST"
echo "=========================================="
echo "Trait: ${TRAIT}${TAG}"
echo "Alignment: $ALI_DIR"
echo "Results: $RESULTS_BASE"
echo "Timestamp: $timestamp"
echo ""

#########################################
# STEP 1: CAAS Post-Processing (Filter)
#########################################
echo "=========================================="
echo "TEST THE WHOLE PIPELINE IN ONE WORKFLOW"
echo "=========================================="
echo "Running complete pipeline: Discovery → Post-processing (filter mode)"
echo ""

RESULTS_BASE_FILTER="${RESULTS_BASE}/filter"
mkdir -p "$RESULTS_BASE_FILTER"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/filter" \
    --prune_data \
    --prune_list "$PRUNE_LIST" \
    --prune_list_secondary "$PRUNE_SECONDARY_LIST" \
    --ct_tool "discovery,resample,bootstrap" \
    --alignment "$ALI_DIR" \
    --my_traits "$TRAIT_FILE" \
    --traitname "$TRAIT" \
    --secondary_trait "$SECONDARY_TRAIT" \
    --branch_trait "$BRANCH_TRAIT" \
    --n_trait "$N_TRAIT" \
    --c_trait "$C_TRAIT" \
    --cycles "$CYCLES" \
    --outdir "$RESULTS_BASE_FILTER" \
    --caas_postproc_mode filter


# echo ""
# echo "✓ Filter mode post-processing completed"
# echo ""

#########################################
# STEP 2: CAAS Post-Processing (Exploratory)
#########################################
echo "=========================================="
echo "TEST THE WHOLE PIPELINE IN ONE WORKFLOW"
echo "=========================================="
echo "Running complete pipeline: Discovery → Post-processing (exploratory mode)"
echo ""

RESULTS_BASE_EXPLORATORY="${RESULTS_BASE}/exploratory"
mkdir -p "$RESULTS_BASE_EXPLORATORY"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/exploratory" \
    --prune_data \
    --prune_list "$PRUNE_LIST" \
    --prune_list_secondary "$PRUNE_SECONDARY_LIST" \
    --ct_tool "discovery,resample,bootstrap" \
    --alignment "$ALI_DIR" \
    --my_traits "$TRAIT_FILE" \
    --traitname "$TRAIT" \
    --secondary_trait "$SECONDARY_TRAIT" \
    --branch_trait "$BRANCH_TRAIT" \
    --n_trait "$N_TRAIT" \
    --c_trait "$C_TRAIT" \
    --cycles "$CYCLES" \
    --outdir "$RESULTS_BASE_EXPLORATORY" \
    --caas_postproc_mode exploratory

#########################################
# Summary
#########################################
echo ""
echo "=========================================="
echo "INTEGRATED PIPELINE TEST COMPLETED"
echo "=========================================="
echo ""
echo "Results directory: $RESULTS_BASE"
echo ""
echo "Pipeline outputs:"
echo "  1. Complete Pipeline (Filter mode):"
echo "     - CAAS Discovery outputs"
echo "     - Resample results"
echo "     - Bootstrap analysis"
echo "     - Post-processing with filter mode"
echo ""
echo "  2. Complete Pipeline (Exploratory mode):"
echo "     - CAAS Discovery outputs"
echo "     - Resample results"
echo "     - Bootstrap analysis"
echo "     - Post-processing with exploratory mode"
echo "     - Parameter sweep results"
echo ""
echo "Key directories to examine:"
echo "  - CAAS results:            ${RESULTS_BASE}/caas/"
echo "  - Resample results:        ${RESULTS_BASE}/resample/"
echo "  - Bootstrap results:       ${RESULTS_BASE}/bootstrap/"
echo "  - Post-processing:         ${RESULTS_BASE}/postprocessing/"
echo ""
echo "Reports and visualizations:"
echo "  - Manhattan plots:         ${RESULTS_BASE}/postprocessing/*/reports/manhattan_plot.png"
echo "  - Characterization:        ${RESULTS_BASE}/postprocessing/*/reports/CAAS_postproc.html"
echo "  - Filter summaries:        ${RESULTS_BASE}/postprocessing/*/discarded_summary.tsv"
echo ""
echo "=========================================="
echo "End-to-end workflow validation complete!"
echo "=========================================="
