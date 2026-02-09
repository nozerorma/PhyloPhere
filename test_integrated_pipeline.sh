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
IS_TOY=true
TRAIT="neoplasia_prevalence"

# Input configuration
if [ "$IS_TOY" = true ]; then
    TAG="_toy"
    ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Ali_toy"
    CYCLES="100"
else
    TAG=""
    ALI_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/Primate_alignments"
    CYCLES="1000000"
fi

# Trait files
TRAIT_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Traitfiles/${TRAIT}/traitfile.tab"
TREE_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/5.Phylogeny/science.abn7829_data_s4.nex.tree"
TRAIT_VALUES="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.5.Bootstrap_traitfiles/${TRAIT}/boot_traitfile.tab"
GENE_ENSEMBL_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/ensembl_genes.output"

# Validate required files
for file in "$TRAIT_FILE" "$TREE_FILE" "$TRAIT_VALUES" "$GENE_ENSEMBL_FILE"; do
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
RESULTS_BASE="Out/integrated_test/${TRAIT}${TAG}/${timestamp}"
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
# STEP 1: CAAS Discovery Pipeline
#########################################
echo "=========================================="
echo "STEP 1: CAAS DISCOVERY PIPELINE"
echo "=========================================="
echo "Running: Discovery + Resample + Bootstrap"
echo ""

OUT_DISCOVERY="${RESULTS_BASE}/caas_analysis"
mkdir -p "$OUT_DISCOVERY"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/discovery" \
    --prune_data \
    --contrast_selection \
    --ct_tool "discovery,resample,bootstrap" \
    --alignment "$ALI_DIR" \
    --caas_config "$TRAIT_FILE" \
    --tree "$TREE_FILE" \
    --traitvalues "$TRAIT_VALUES" \
    --cycles "$CYCLES" \
    --outdir "$OUT_DISCOVERY"

echo ""
echo "✓ CAAS discovery completed"
echo ""

# Locate output files from discovery
DISCOVERY_FILE=$(find "$OUT_DISCOVERY" -name "discovery.tab" | head -1)
BACKGROUND_FILE=$(find "$OUT_DISCOVERY" -name "background_genes.output" | head -1)

if [ ! -f "$DISCOVERY_FILE" ]; then
    echo "Error: Discovery file not found in $OUT_DISCOVERY"
    echo "Expected: discovery.tab"
    exit 1
fi

if [ ! -f "$BACKGROUND_FILE" ]; then
    echo "Error: Background genes file not found in $OUT_DISCOVERY"
    echo "Expected: background_genes.output"
    exit 1
fi

echo "Discovery outputs located:"
echo "  Discovery: $DISCOVERY_FILE"
echo "  Background: $BACKGROUND_FILE"
echo "  Gene ensembl: $GENE_ENSEMBL_FILE"
echo ""

#########################################
# STEP 2: CAAS Post-Processing (Exploratory)
#########################################
echo "=========================================="
echo "STEP 2: CAAS POST-PROCESSING (EXPLORATORY)"
echo "=========================================="
echo "Running parameter sweep: minlen=2,3,4 × maxcaas=0.6,0.7,0.8"
echo ""

OUT_POSTPROC_EXPLORATORY="${RESULTS_BASE}/postprocessing/exploratory"
mkdir -p "$OUT_POSTPROC_EXPLORATORY"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/postproc_exploratory" \
    --caas_postproc \
    --discovery_input "$DISCOVERY_FILE" \
    --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
    --background_input "$BACKGROUND_FILE" \
    --caas_postproc_mode exploratory \
    --minlen_values "2,3,4" \
    --maxcaas_values "0.6,0.7,0.8" \
    --gene_filter_mode both \
    --extreme_threshold 0.99 \
    --iqr_multiplier 3.0 \
    --generate_reports true \
    --generate_manhattan true \
    --postproc_outdir "$OUT_POSTPROC_EXPLORATORY" \
    --outdir "$OUT_POSTPROC_EXPLORATORY"

echo ""
echo "✓ Exploratory post-processing completed"
echo ""

#########################################
# STEP 3: CAAS Post-Processing (Filter)
#########################################
echo "=========================================="
echo "STEP 3: CAAS POST-PROCESSING (FILTER)"
echo "=========================================="
echo "Applying best-performing filter: minlen=2, maxcaas=0.6"
echo ""

OUT_POSTPROC_FILTER="${RESULTS_BASE}/postprocessing/filter"
mkdir -p "$OUT_POSTPROC_FILTER"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/postproc_filter" \
    --caas_postproc \
    --discovery_input "$DISCOVERY_FILE" \
    --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
    --background_input "$BACKGROUND_FILE" \
    --caas_postproc_mode filter \
    --filter_minlen 2 \
    --filter_maxcaas 0.6 \
    --gene_filter_mode both \
    --extreme_threshold 0.99 \
    --iqr_multiplier 3.0 \
    --generate_reports true \
    --generate_manhattan true \
    --manhattan_min_density 0.008 \
    --postproc_outdir "$OUT_POSTPROC_FILTER" \
    --outdir "$OUT_POSTPROC_FILTER"

echo ""
echo "✓ Filter mode post-processing completed"
echo ""

#########################################
# STEP 4: Test Different Gene Filter Modes
#########################################
echo "=========================================="
echo "STEP 4: GENE FILTER MODE COMPARISON"
echo "=========================================="
echo "Testing: extreme, dubious, both, none"
echo ""

for filter_mode in extreme dubious both none; do
    echo "→ Testing gene_filter_mode=${filter_mode}"
    
    OUT_FILTER_MODE="${RESULTS_BASE}/postprocessing/gene_filter_${filter_mode}"
    mkdir -p "$OUT_FILTER_MODE"
    
    nextflow run main.nf \
        -with-tower \
        -profile local \
        -w "$WORK_DIR/gene_filter_${filter_mode}" \
        --caas_postproc \
        --discovery_input "$DISCOVERY_FILE" \
        --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
        --background_input "$BACKGROUND_FILE" \
        --caas_postproc_mode filter \
        --filter_minlen 2 \
        --filter_maxcaas 0.6 \
        --gene_filter_mode "$filter_mode" \
        --extreme_threshold 0.99 \
        --iqr_multiplier 3.0 \
        --generate_reports true \
        --generate_manhattan true \
        --postproc_outdir "$OUT_FILTER_MODE" \
        --outdir "$OUT_FILTER_MODE"
    
    echo "  ✓ ${filter_mode} completed"
    echo ""
done

#########################################
# Summary
#########################################
echo "=========================================="
echo "INTEGRATED PIPELINE TEST COMPLETED"
echo "=========================================="
echo ""
echo "Results directory: $RESULTS_BASE"
echo ""
echo "Pipeline outputs:"
echo "  1. CAAS Analysis:           $OUT_DISCOVERY"
echo "     - Discovery file:        $(basename "$DISCOVERY_FILE")"
echo "     - Background file:       $(basename "$BACKGROUND_FILE")"
echo ""
echo "  2. Post-processing (Exploratory):"
echo "     - Directory:             $OUT_POSTPROC_EXPLORATORY"
echo "     - Parameter sweep:       9 combinations tested"
echo "     - Reports:               */reports/*.html"
echo ""
echo "  3. Post-processing (Filter):"
echo "     - Directory:             $OUT_POSTPROC_FILTER"
echo "     - Best filter:           minlen=2, maxcaas=0.6"
echo "     - Filtered output:       filter_*/discovery.filtered.*.tsv"
echo ""
echo "  4. Gene Filter Comparison:"
echo "     - Extreme only:          ${RESULTS_BASE}/postprocessing/gene_filter_extreme"
echo "     - Dubious only:          ${RESULTS_BASE}/postprocessing/gene_filter_dubious"
echo "     - Both filters:          ${RESULTS_BASE}/postprocessing/gene_filter_both"
echo "     - No filtering:          ${RESULTS_BASE}/postprocessing/gene_filter_none"
echo ""
echo "Key files to examine:"
echo "  - Filter summaries:        */discarded_summary.tsv"
echo "  - Gene filtering:          */gene_filtering/removed_genes_summary.tsv"
echo "  - Manhattan plots:         */reports/manhattan_plot.png"
echo "  - Characterization:        */reports/CAAS_postproc.html"
echo ""
echo "=========================================="
echo "End-to-end workflow validation complete!"
echo "=========================================="
