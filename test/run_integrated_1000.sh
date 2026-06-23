#!/usr/bin/env bash

set -Eeuo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="${REPO_DIR}/test"
INPUTS_DIR="${TEST_DIR}/inputs"
ALIGNMENT_DIR="${REPO_DIR}/test/inputs/alignments/Ali_1000"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${TEST_DIR}/runs/integrated_1000/${TIMESTAMP}"
WORKDIR="${OUTDIR}/work"

TRAIT_FILE="${INPUTS_DIR}/traitfiles/traitfile.tab"
BOOT_TRAIT_FILE="${INPUTS_DIR}/traitfiles/boot_traitfile.tab"
TREE_FILE="${INPUTS_DIR}/tree/pruned_tree_file.nwk"
ASR_CACHE_DIR="${REPO_DIR}/test/inputs/asr"

require_path() {
    local path="$1"
    if [[ ! -e "$path" ]]; then
        echo "Error: required path not found: $path" >&2
        exit 1
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# ENVIRONMENT ACTIVATION
# ─────────────────────────────────────────────────────────────────────────────
echo "Activating environment: phylophere"
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
elif command -v conda &> /dev/null; then
    conda activate phylophere || echo "Warning: Could not activate conda environment"
else
    echo "Warning: Neither micromamba nor conda found — proceeding without activation"
fi

require_path "$ALIGNMENT_DIR"
require_path "$TRAIT_FILE"
require_path "$BOOT_TRAIT_FILE"
require_path "$TREE_FILE"
require_path "$ASR_CACHE_DIR"

mkdir -p "$OUTDIR" "$WORKDIR"

echo "=========================================="
echo "PHYLOPHERE INTEGRATED 1000 GENES TESTER"
echo "=========================================="
echo "Mode              : full integrated CT 1000 run"
echo "Alignments        : $ALIGNMENT_DIR"
echo "Trait file        : $TRAIT_FILE"
echo "Boot trait file   : $BOOT_TRAIT_FILE"
echo "Tree file         : $TREE_FILE"
echo "ASR cache         : $ASR_CACHE_DIR"
echo "Output dir        : $OUTDIR"
echo "Work dir          : $WORKDIR"
echo "=========================================="

NF_FLAGS=(
    -profile local
    --ali_format fasta
    --reporting false
    --ct_tool "discovery,resample,bootstrap"
    --ct_signification
    --ct_disambiguation
    --ct_postproc
    --ct_accumulation
    --traitname "neoplasia_prevalence"
    --scoring
    --gene_ensembl_file "${REPO_DIR}/test/inputs/alignments/ensembl_genes.output"
    --tax_id "${REPO_DIR}/test/inputs/phylogeny/taxid_species_family_primates.tsv"
    --vep_map_dir "${REPO_DIR}/test/inputs/map_1000/"
    --alignment "$ALIGNMENT_DIR"
    --caas_config "$TRAIT_FILE"
    --traitvalues "$BOOT_TRAIT_FILE"
    --tree "$TREE_FILE"
    --ct_disambig_asr_mode precomputed
    --ct_disambig_asr_cache_dir "$ASR_CACHE_DIR"
    --ct_disambig_posterior_threshold 0.1
    --caas_postproc_mode filter
    --gene_filter_mode both
    --generate_reports true
    --asr_robustness true
    --scoring_rer_input "${REPO_DIR}/test/inputs/rer/RERConverge/RER_Results/rerconverge_summary_neoplasia_prevalence.tsv"
    --vep
    --outdir "$OUTDIR"
)

nextflow run "${REPO_DIR}/main.nf" \
    "${NF_FLAGS[@]}" \
    -w "$WORKDIR" \
    "$@"

echo ""
echo "Integrated 1000 genes test finished."
echo "Results: $OUTDIR"
