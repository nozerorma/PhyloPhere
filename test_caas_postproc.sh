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
# PHYLOPHERE: CAAS Post-Processing Test Script
#
# Tests CAAS post-processing workflow in both exploratory and filter modes
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: test_caas_postproc.sh
#

set -Eeuo pipefail

# Activate micromamba environment
echo "Activating micromamba environment: caas_global_cancer"
# Initialize micromamba for bash
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
# Try to find micromamba
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate caas_global_cancer
else
    echo "Warning: micromamba not found in PATH, attempting manual conda activation"
    # Try conda as fallback
    if command -v conda &> /dev/null; then
        conda activate caas_global_cancer || echo "Warning: Could not activate conda environment"
    else
        echo "Warning: Neither micromamba nor conda found, proceeding without environment activation"
    fi
fi

timestamp=$(date +%Y%m%d_%H%M%S)

# Test configuration
TEST_TAG="postproc_test"
TRAIT="neoplasia_prevalence_toy"

# Input files (using existing toy dataset from PhyloPhere Out directory)
DISCOVERY_FILE="Out/integrated_test/neoplasia_prevalence_toy/20260210_095707/caas_analysis/caastools/discovery.tab"
BOOTSTRAP_FILE="Out/integrated_test/neoplasia_prevalence_toy/20260210_095707/caas_analysis/caastools/bootstrap.tab"
GENE_ENSEMBL_FILE="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data/2.Alignments/ensembl_genes.output"
BACKGROUND_FILE="Out/integrated_test/neoplasia_prevalence_toy/20260210_095707/caas_analysis/caastools/background_genes.output"

# Validate input files exist
if [ ! -f "$DISCOVERY_FILE" ]; then
    echo "Error: Discovery file not found: $DISCOVERY_FILE"
    exit 1
fi

if [ ! -f "$BOOTSTRAP_FILE" ]; then
    echo "Error: Bootstrap file not found: $BOOTSTRAP_FILE"
    exit 1
fi

if [ ! -f "$GENE_ENSEMBL_FILE" ]; then
    echo "Error: Gene ensembl file not found: $GENE_ENSEMBL_FILE"
    exit 1
fi

if [ ! -f "$BACKGROUND_FILE" ]; then
    echo "Error: Background file not found: $BACKGROUND_FILE"
    echo "Background files are required for CAAS post-processing (gene filtering cleanup)"
    exit 1
fi

echo "=========================================="
echo "CAAS Post-Processing Test"
echo "=========================================="
echo "Trait: $TRAIT"
echo "Discovery: $DISCOVERY_FILE"
echo "Bootstrap: $BOOTSTRAP_FILE"
echo "Gene ensembl: $GENE_ENSEMBL_FILE"
echo "Background: $BACKGROUND_FILE"
echo "Timestamp: $timestamp"
echo ""

# Create output directory structure
RESULTS_BASE="Out/caas_postproc_test/${TRAIT}/${timestamp}"
mkdir -p "$RESULTS_BASE"

WORK_DIR="${RESULTS_BASE}/work"
mkdir -p "$WORK_DIR"

# For Nextflow Tower/Apptainer
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_SINGULARITY_HOME_MOUNT=true

#########################################
# TEST 1: Exploratory Mode
#########################################
echo "=========================================="
echo "TEST 1: EXPLORATORY MODE"
echo "=========================================="
echo "Testing parameter sweep (2,3,4 × 0.6,0.7,0.8)"
echo ""

OUT_EXPLORATORY="${RESULTS_BASE}/exploratory/postproc"
OUT_SIGNIFICATION="${RESULTS_BASE}/exploratory/caas_signification"
mkdir -p "$OUT_EXPLORATORY"
mkdir -p "$OUT_SIGNIFICATION"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/exploratory" \
    --caas_postproc \
    --caas_signification \
    --discovery_input "$DISCOVERY_FILE" \
    --bootstrap_input "$BOOTSTRAP_FILE" \
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
    --postproc_outdir "$OUT_EXPLORATORY" \
    --signification_outdir "$OUT_SIGNIFICATION"

echo ""
echo "✓ Exploratory mode completed"
echo "  Results: $OUT_EXPLORATORY"
echo ""

#########################################
# TEST 2: Filter Mode
#########################################
echo "=========================================="
echo "TEST 2: FILTER MODE"
echo "=========================================="
echo "Testing single filter (minlen=3, maxcaas=0.7)"
echo ""

OUT_FILTER="${RESULTS_BASE}/filter/postproc"
mkdir -p "$OUT_FILTER"
OUT_SIGNIFICATION="${RESULTS_BASE}/filter/caas_signification"
mkdir -p "$OUT_SIGNIFICATION"

nextflow run main.nf \
    -with-tower \
    -profile local \
    -w "$WORK_DIR/filter" \
    --caas_postproc \
    --caas_signification \
    --discovery_input "$DISCOVERY_FILE" \
    --bootstrap_input "$BOOTSTRAP_FILE" \
    --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
    --background_input "$BACKGROUND_FILE" \
    --caas_postproc_mode filter \
    --filter_minlen 3 \
    --filter_maxcaas 0.7 \
    --gene_filter_mode both \
    --extreme_threshold 0.99 \
    --iqr_multiplier 3.0 \
    --generate_reports true \
    --generate_manhattan true \
    --postproc_outdir "$OUT_FILTER" \
    --signification_outdir "$OUT_SIGNIFICATION"

echo ""
echo "✓ Filter mode completed"
echo "  Results: $OUT_FILTER"
echo ""

# #########################################
# # TEST 3: Gene Filtering Modes
# #########################################
# echo "=========================================="
# echo "TEST 3: GENE FILTERING MODES"
# echo "=========================================="
# echo "Testing different gene filter modes"
# echo ""

# # Test with only extreme filtering
# OUT_EXTREME="${RESULTS_BASE}/gene_filter_extreme"
# mkdir -p "$OUT_EXTREME"
# OUT_SIGNIFICATION="${RESULTS_BASE}/gene_filter_extreme/caas_signification"
# mkdir -p "$OUT_SIGNIFICATION"

# echo "→ Testing gene_filter_mode=extreme"
# nextflow run main.nf \
#     -with-tower \
#     -profile local \
#     -w "$WORK_DIR/gene_extreme" \
#     --caas_postproc \
#     --caas_signification \
#     --discovery_input "$DISCOVERY_FILE" \
#     --bootstrap_input "$BOOTSTRAP_FILE" \
#     --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
#     --background_input "$BACKGROUND_FILE" \
#     --caas_postproc_mode filter \
#     --filter_minlen 3 \
#     --filter_maxcaas 0.7 \
#     --gene_filter_mode extreme \
#     --extreme_threshold 0.95 \
#     --generate_reports true \
#     --postproc_outdir "$OUT_EXTREME" \
#     --signification_outdir "$OUT_SIGNIFICATION"

# echo "  ✓ Extreme-only filtering completed"
# echo ""

# # Test with only dubious filtering
# OUT_DUBIOUS="${RESULTS_BASE}/gene_filter_dubious"
# mkdir -p "$OUT_DUBIOUS"
# OUT_SIGNIFICATION="${RESULTS_BASE}/gene_filter_dubious/caas_signification"
# mkdir -p "$OUT_SIGNIFICATION"

# echo "→ Testing gene_filter_mode=dubious"
# nextflow run main.nf \
#     -with-tower \
#     -profile local \
#     -w "$WORK_DIR/gene_dubious" \
#     --caas_postproc \
#     --caas_signification \
#     --discovery_input "$DISCOVERY_FILE" \
#     --bootstrap_input "$BOOTSTRAP_FILE" \
#     --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
#     --background_input "$BACKGROUND_FILE" \
#     --caas_postproc_mode filter \
#     --filter_minlen 3 \
#     --filter_maxcaas 0.7 \
#     --gene_filter_mode dubious \
#     --iqr_multiplier 4.0 \
#     --generate_reports true \
#     --postproc_outdir "$OUT_DUBIOUS" \
#     --signification_outdir "$OUT_SIGNIFICATION"

# echo "  ✓ Dubious-only filtering completed"
# echo ""

# # Test with no gene filtering
# OUT_NONE="${RESULTS_BASE}/gene_filter_none"
# mkdir -p "$OUT_NONE"
# OUT_SIGNIFICATION="${RESULTS_BASE}/gene_filter_none/caas_signification"
# mkdir -p "$OUT_SIGNIFICATION"

# echo "→ Testing gene_filter_mode=none"
# nextflow run main.nf \
#     -with-tower \
#     -profile local \
#     -w "$WORK_DIR/gene_none" \
#     --caas_postproc \
#     --caas_signification \
#     --discovery_input "$DISCOVERY_FILE" \
#     --bootstrap_input "$BOOTSTRAP_FILE" \
#     --gene_ensembl_file "$GENE_ENSEMBL_FILE" \
#     --background_input "$BACKGROUND_FILE" \
#     --caas_postproc_mode filter \
#     --filter_minlen 3 \
#     --filter_maxcaas 0.7 \
#     --gene_filter_mode none \
#     --generate_reports true \
#     --postproc_outdir "$OUT_NONE" \
#     --signification_outdir "$OUT_SIGNIFICATION"

# echo "  ✓ No gene filtering completed"
# echo ""

#########################################
# Summary
#########################################
echo "=========================================="
echo "ALL TESTS COMPLETED"
echo "=========================================="
echo ""
echo "Results summary:"
echo "  Base directory: $RESULTS_BASE"
echo ""
echo "  1. Exploratory mode:      $OUT_EXPLORATORY"
echo "  2. Filter mode:           $OUT_FILTER"
echo "  3. Extreme-only filter:   $OUT_EXTREME"
echo "  4. Dubious-only filter:   $OUT_DUBIOUS"
echo "  5. No gene filter:        $OUT_NONE"
echo ""
echo "Check the following for outputs:"
echo "  - Filtered CAAS files:    */filter_*/*.filtered.*.tsv"
echo "  - Filter summaries:       */discarded_summary.tsv"
echo "  - Gene filtering results: */gene_filtering/"
echo "  - Reports:                */reports/*.html"
echo "  - Manhattan plots:        */reports/manhattan_plot.png"
echo ""
echo "=========================================="
echo "Test script completed successfully!"
echo "=========================================="
