#!/bin/bash

# Quick validation script for CAAS Post-Processing workflow integration
# This script checks that all required files are in place and properly configured

set -euo pipefail

BASE_DIR="/home/miguel/IBE-UPF/PhD/PhyloPhere"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  CAAS Post-Processing Workflow - Integration Validation       ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check main workflow file
echo "✓ Checking main workflow file..."
if [ -f "$BASE_DIR/workflows/caas_postproc.nf" ]; then
    echo "  ✓ workflows/caas_postproc.nf exists"
else
    echo "  ✗ workflows/caas_postproc.nf NOT FOUND"
    exit 1
fi

# Check subworkflow files
echo "✓ Checking subworkflow files..."
for file in \
    "subworkflows/CAAS_POSTPROC/ctpp_clustfilter.nf" \
    "subworkflows/CAAS_POSTPROC/ctpp_characterization.nf" \
    "subworkflows/CAAS_POSTPROC/README.md"; do
    if [ -f "$BASE_DIR/$file" ]; then
        echo "  ✓ $file exists"
    else
        echo "  ✗ $file NOT FOUND"
        exit 1
    fi
done

# Check local supporting files
echo "✓ Checking supporting files..."
for file in \
    "subworkflows/CAAS_POSTPROC/local/filter_caas_clusters-param.py" \
    "subworkflows/CAAS_POSTPROC/local/CAAS_postproc.Rmd"; do
    if [ -f "$BASE_DIR/$file" ]; then
        echo "  ✓ $file exists"
    else
        echo "  ✗ $file NOT FOUND"
        exit 1
    fi
done

# Check configuration file
echo "✓ Checking configuration file..."
if [ -f "$BASE_DIR/conf/caas_postproc.config" ]; then
    echo "  ✓ conf/caas_postproc.config exists"
else
    echo "  ✗ conf/caas_postproc.config NOT FOUND"
    exit 1
fi

# Check main.nf includes
echo "✓ Verifying main.nf integration..."
if grep -q "include {CAAS_POSTPROC}" "$BASE_DIR/main.nf"; then
    echo "  ✓ CAAS_POSTPROC workflow included in main.nf"
else
    echo "  ✗ CAAS_POSTPROC workflow NOT included in main.nf"
    exit 1
fi

if grep -q "params.caas_postproc" "$BASE_DIR/main.nf"; then
    echo "  ✓ caas_postproc parameter check in main.nf"
else
    echo "  ✗ caas_postproc parameter check NOT found in main.nf"
    exit 1
fi

# Check nextflow.config includes
echo "✓ Verifying nextflow.config integration..."
if grep -q "includeConfig 'conf/caas_postproc.config'" "$BASE_DIR/nextflow.config"; then
    echo "  ✓ caas_postproc.config included in nextflow.config"
else
    echo "  ✗ caas_postproc.config NOT included in nextflow.config"
    exit 1
fi

# Check common.config
echo "✓ Verifying common.config integration..."
if grep -q "caas_postproc" "$BASE_DIR/conf/common.config"; then
    echo "  ✓ caas_postproc parameter defined in common.config"
else
    echo "  ✗ caas_postproc parameter NOT defined in common.config"
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  ✅ All validation checks passed!                             ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Integration Summary:"
echo "  • Main orchestrator: workflows/caas_postproc.nf"
echo "  • Cluster filtering: subworkflows/CAAS_POSTPROC/ctpp_clustfilter.nf"
echo "  • Characterization: subworkflows/CAAS_POSTPROC/ctpp_characterization.nf"
echo "  • Configuration: conf/caas_postproc.config"
echo ""
echo "Usage Examples:"
echo "  1. Exploratory mode (parameter sweep):"
echo "     nextflow run main.nf --caas_postproc --discovery_input data/discovery.tab --gene_length_file data/genes.tsv"
echo ""
echo "  2. Filter mode (single parameters):"
echo "     nextflow run main.nf --caas_postproc --caas_postproc_mode filter --filter_minlen 3 --filter_maxcaas 0.7 ..."
echo ""
echo "  3. Without reports:"
echo "     nextflow run main.nf --caas_postproc --discovery_input data/discovery.tab --generate_reports false"
echo ""
echo "For detailed documentation, see: subworkflows/CAAS_POSTPROC/README.md"
