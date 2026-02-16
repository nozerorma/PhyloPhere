#!/usr/bin/env bash
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
#
# run_pipeline.sh - Simple orchestrator for single-gene CAAS analysis
# Follows the pattern of run_aggregate.sh for simplicity and maintainability.
#
# Usage:
#   ./run_pipeline.sh -g GENE [GENE...] [OPTIONS]
#
# Options:
#   -g, --genes GENE [GENE...]    Gene names to analyze (required)
#   -a, --alignment-dir DIR       Alignment directory (default: relative path)
#   -t, --tree FILE              Tree file (default: relative path)
#   -o, --output DIR            Output directory (default: timestamped)
#   --taxid-mapping FILE        TaxID mapping file
#   --caas-metadata FILE        CAAS metadata file
#   --trait-file FILE           Trait file
#   --caas-positions POS...     CAAS positions (space-separated)
#   --asr-mode MODE            ASR mode: compute | precomputed (default: precomputed)
#   --with-asr                  Legacy alias, forces --asr-mode compute
#   --asr-model MODEL           ASR model (default: lg)
#   --asr-cache-dir DIR         Directory for ASR cache files
#   --ensembl-genes-file FILE   TSV/CSV with gene column to enforce exact alignment matching
#   --convergence-mode MODE     Convergence mode: 'focal_clade' (default) or 'mrca'
#   --run-diagnostics           Enable diagnostics: node-level posteriors (TSV) and tip details (JSON) in diagnostics/
#   --focus-positions           Highlight three sampled positions in plots
#   --debug-positions           Enable full debug outputs for first three analyzed sites
#   --allow-low-confidence      Continue even when ASR confidence is low
#   --include-non-significant   Include non-significant CAAS in reports
#   --threads N                 Threads for codeml (default: 1)
#   -v, --verbose               Verbose logging
#   -h, --help                  Show this help
#
# Examples:
#   # Analyze NUTM2A with default settings
#   ./run_pipeline.sh -g NUTM2A
#
#   # Multiple genes with ASR
#   ./run_pipeline.sh -g BRCA1 TP53 --with-asr
#
#   # Custom positions and output
#   ./run_pipeline.sh -g NUTM2A --caas-positions 85 176 206 -o my_analysis

set -Eeuo pipefail

# -----------------------------
# Defaults (relative to script location)
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data"
RESULTS_DIR="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out"
ALIGNMENT_DIR="${DATA_DIR}/2.Alignments/Primate_alignments"
TREE_FILE="${DATA_DIR}/5.Phylogeny/science.abn7829_data_s4.nex.tree"
TAXID_MAPPING="${DATA_DIR}/5.Phylogeny/taxid_species_family.tsv"
TRAIT_FILE="${RESULTS_DIR}/2.CAAS/1.Discovery/1.Traitfiles/malignant_prevalence/malignant_prevalence_paired_species.tab"
PIPELINE_SCRIPT="${SCRIPT_DIR}/pipeline_main.py"

# Optional defaults
ASR_MODEL="${ASR_MODEL:-lg}"
OUTPUT_DIR="${SCRIPT_DIR}/results/single_gene_pipeline"
GENES="TP53"
CAAS_POSITIONS=""
POSTERIOR_THRESHOLD="0.0"
CAAS_METADATA="${RESULTS_DIR}/2.CAAS/3.5.CAAS_exploration/4.Commonalities/not-significant/malignant_prevalence/meta_caas.output"
WITH_ASR=""
ASR_MODE="precomputed"
ASR_CACHE_DIR="${SCRIPT_DIR}/results/asr_cache"
VERBOSE=""
CONVERGENCE_MODE="mrca"
RUN_DIAGNOSTICS=""
FOCUS_POSITIONS="--focus-positions"
DEBUG_POSITIONS="--debug-positions"
ALLOW_LOW_CONF="--allow-low-confidence"
INCLUDE_NON_SIGNIFICANT=""
THREADS="4"
ENSEMBL_GENES_FILE=""

# -----------------------------
# Helpers
# -----------------------------
usage() {
  sed -n '3,25p' "$0"
  exit 1
}

log() { printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*"; }

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    log "ERROR: $cmd not found"
    exit 3
  fi
}

# Select Python runner (prefer micromamba env if available)
PY_RUNNER="python3"
MAMBA_ENV="${MAMBA_ENV:-caas_global_cancer}"
if command -v micromamba >/dev/null 2>&1; then
  PY_RUNNER="micromamba run -n ${MAMBA_ENV} python"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/..:${PYTHONPATH:-}"

ASR_MODE=""

# -----------------------------
# Parse CLI
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genes) GENES="$2"; shift 2;;
    -a|--alignment-dir) ALIGNMENT_DIR="$2"; shift 2;;
    -t|--tree) TREE_FILE="$2"; shift 2;;
    -o|--output) OUTPUT_DIR="$2"; shift 2;;
    --taxid-mapping) TAXID_MAPPING="$2"; shift 2;;
    --caas-metadata) CAAS_METADATA="$2"; shift 2;;
    --trait-file) TRAIT_FILE="$2"; shift 2;;
    --caas-positions)
      shift
      while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
        CAAS_POSITIONS="$CAAS_POSITIONS $1"
        shift
      done
      ;;
    --with-asr) ASR_MODE="compute"; shift;;
    --asr-mode) ASR_MODE="$2"; shift 2;;
    --asr-cache-dir) ASR_CACHE_DIR="$2"; shift 2;;
    --asr-model) ASR_MODEL="$2"; shift 2;;
    --ensembl-genes-file) ENSEMBL_GENES_FILE="$2"; shift 2;;
    --convergence-mode) CONVERGENCE_MODE="$2"; shift 2;;
    --posterior-threshold) POSTERIOR_THRESHOLD="$2"; shift 2;;
    --run-diagnostics) RUN_DIAGNOSTICS="--run-diagnostics"; shift;;
    --focus-positions) FOCUS_POSITIONS="--focus-positions"; shift;;
    --debug-positions) DEBUG_POSITIONS="--debug-positions"; shift;;
    --allow-low-confidence) ALLOW_LOW_CONF="--allow-low-confidence"; shift;;
    --include-non-significant) INCLUDE_NON_SIGNIFICANT="--include-non-significant"; shift;;
    --threads) THREADS="$2"; shift 2;;
    -v | --verbose) VERBOSE="--verbose"; shift;;
    -h | --help) usage;;
    *) log "Unknown option: $1"; usage;;
  esac
done

# -----------------------------
# Validation
# -----------------------------
require_cmd ${PY_RUNNER%% *}

if [[ -z "$GENES" ]]; then
  log "ERROR: --genes is required"
  usage
fi

if [[ ! -d "$ALIGNMENT_DIR" ]]; then
  log "ERROR: Alignment directory not found: $ALIGNMENT_DIR"
  exit 1
fi

if [[ ! -f "$TREE_FILE" ]]; then
  log "ERROR: Tree file not found: $TREE_FILE"
  exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="./pipeline_output_$(date +%Y%m%d_%H%M%S)"
fi

# Create output dir
mkdir -p "$OUTPUT_DIR"

# -----------------------------
# Navigation and execution
# -----------------------------
# Navigation and execution (run from parent directory for relative imports)
pushd "$SCRIPT_DIR/.." > /dev/null

log "Starting CAAS pipeline for genes: $GENES"
log "Alignment dir: $ALIGNMENT_DIR"
log "Tree: $TREE_FILE"
log "Output: $OUTPUT_DIR"

# Build command
cmd="${PY_RUNNER} -m single_gene_pipeline.pipeline_main"
cmd="$cmd --genes $GENES"
cmd="$cmd --alignment_dir \"$ALIGNMENT_DIR\""
cmd="$cmd --tree \"$TREE_FILE\""
cmd="$cmd --output_dir \"$OUTPUT_DIR\""
cmd="$cmd --posterior-threshold $POSTERIOR_THRESHOLD"
cmd="$cmd --convergence-mode \"$CONVERGENCE_MODE\""
cmd="$cmd --threads $THREADS"

[[ -n "$CAAS_POSITIONS" ]] && cmd="$cmd --caas_positions $CAAS_POSITIONS"
[[ -n "$CAAS_METADATA" ]] && cmd="$cmd --caas_metadata \"$CAAS_METADATA\""
[[ -n "$TAXID_MAPPING" ]] && cmd="$cmd --taxid_mapping \"$TAXID_MAPPING\""
[[ -n "$TRAIT_FILE" ]] && cmd="$cmd --trait_file \"$TRAIT_FILE\""
[[ -n "$ASR_MODE" ]] && cmd="$cmd --asr-mode $ASR_MODE"
[[ -n "$ASR_CACHE_DIR" ]] && cmd="$cmd --asr-cache-dir \"$ASR_CACHE_DIR\""
[[ -n "$ENSEMBL_GENES_FILE" ]] && cmd="$cmd --ensembl-genes-file \"$ENSEMBL_GENES_FILE\""
[[ -n "$VERBOSE" ]] && cmd="$cmd $VERBOSE"
cmd="$cmd --asr_model \"$ASR_MODEL\""
[[ -n "$RUN_DIAGNOSTICS" ]] && cmd="$cmd $RUN_DIAGNOSTICS"
[[ -n "$FOCUS_POSITIONS" ]] && cmd="$cmd $FOCUS_POSITIONS"
[[ -n "$DEBUG_POSITIONS" ]] && cmd="$cmd $DEBUG_POSITIONS"
[[ -n "$ALLOW_LOW_CONF" ]] && cmd="$cmd $ALLOW_LOW_CONF"
[[ -n "$INCLUDE_NON_SIGNIFICANT" ]] && cmd="$cmd $INCLUDE_NON_SIGNIFICANT"

# Execute
log "Command: $cmd"
eval "$cmd"

popd > /dev/null
