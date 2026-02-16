#!/usr/bin/env bash
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
#
# run_disambiguation.sh - Orchestrator for CAAS disambiguation across all genes
# Processes all genes found in metadata file and produces consolidated CSVs.
# Toy mode now mirrors 6.CAAS_significance_test/run_aggregate.sh: selecting toy
# alignments when requested, auto-tagging outputs with _toy, and falling back to
# real alignments if toy data are absent.
#
# Usage:
#   ./run_disambiguation.sh [OPTIONS]
#
# Options (key ones):
#   -a, --alignment-dir DIR      Alignment directory (default: Primate_alignments)
#   --alignments-toy DIR         Toy alignment directory (default: Ali_toy)
#   -t, --tree FILE              Tree file (default: science.abn7829_data_s4.nex.tree)
#   -o, --output-dir DIR         Output directory (default: results/<trait>[_toy]_TIMESTAMP)
#   --toy true|false             Use TOY mode (default: true)
#   --tag SUFFIX                 Custom tag for output files (auto: _toy when --toy true)
#   --taxid-mapping FILE         TaxID mapping file (default: taxid_species_family.tsv)
#   --caas-metadata FILE         CAAS metadata file (default: meta_caas.output for trait)
#   --trait-file FILE            Trait file (default: <trait>_paired_species.tab)
#   --ensembl-genes-file FILE    Ensembl genes file (default: ensembl_genes.output)
#   --asr-mode MODE              ASR mode: compute | precomputed (default: precomputed)
#   --asr-model MODEL            ASR model (default: lg)
#   --asr-cache-dir DIR          Directory for ASR cache files
#   --convergence-mode MODE      Convergence mode: 'focal_clade' or 'mrca' (default: mrca)
#   --posterior-threshold VAL    Posterior probability threshold (default: 0.0)
#   --run-diagnostics            Enable diagnostics output: node-level posteriors and tip residues to diagnostics/
#   --skip-gene-lists            Skip per-pattern gene list export (gene_reports/)
#   --allow-low-confidence       Continue even when ASR confidence is low
#   --include-non-significant    Include non-significant CAAS in outputs
#   --threads N                  Threads for codeml per gene (default: 1)
#   --workers N                  Number of parallel gene workers (default: cpu_count-2)
#   -v, --verbose                Verbose logging
#   -h, --help                   Show this help

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TRAIT="${TRAIT:-malignant_prevalence}"
TOY="${TOY:-false}"

if $TOY; then
  TAG="_toy"
else
  TAG=""
fi

DATA_ROOT="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Data"
RESULTS_ROOT="/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out"

ALIGN_ROOT="${ALIGN_ROOT:-${DATA_ROOT}/2.Alignments}"
ALIGNMENT_DIR_REAL_DEFAULT="${ALIGNMENT_DIR_REAL_DEFAULT:-${ALIGN_ROOT}/Primate_alignments}"
ALIGNMENT_DIR_TOY_DEFAULT="${ALIGNMENT_DIR_TOY_DEFAULT:-${ALIGN_ROOT}/Ali_toy}"
ALIGNMENT_DIR_REAL="${ALIGNMENT_DIR_REAL:-${ALIGNMENT_DIR_REAL_DEFAULT}}"
ALIGNMENT_DIR_TOY="${ALIGNMENT_DIR_TOY:-${ALIGNMENT_DIR_TOY_DEFAULT}}"

PHYLO_DIR="${PHYLO_DIR:-${DATA_ROOT}/5.Phylogeny}"
TREE_FILE="${TREE_FILE:-${PHYLO_DIR}/science.abn7829_data_s4.nex.tree}"
TAXID_MAPPING="${TAXID_MAPPING:-${PHYLO_DIR}/taxid_species_family.tsv}"

CAAS_METADATA="${CAAS_METADATA:-${RESULTS_ROOT}/2.CAAS/3.5.CAAS_exploration/4.Commonalities/not-significant/${TRAIT}/meta_caas.output}"
TRAIT_FILE="${TRAIT_FILE:-${RESULTS_ROOT}/2.CAAS/1.Discovery/1.Traitfiles/${TRAIT}/${TRAIT}_paired_species.tab}"
ENSEMBL_GENES_FILE="${ENSEMBL_GENES_FILE:-${ALIGN_ROOT}/ensembl_genes.output}"

ASR_MODEL="${ASR_MODEL:-lg}"
OUTPUT_DIR="${OUTPUT_DIR:-/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92//3.CAAS_disambiguation/${TRAIT}${TAG}}"
POSTERIOR_THRESHOLD="0.7"
ASR_MODE="${ASR_MODE:-precomputed}"
ASR_CACHE_DIR="${ASR_CACHE_DIR:-/media/miguel/adfbf391-5867-414b-8af7-bceb102e6e92/asr}"
CONVERGENCE_MODE="focal_clade"
THREADS="4"
WORKERS="4"
RUN_DIAGNOSTICS="--run-diagnostics"
SKIP_GENE_LISTS=""
INCLUDE_NON_SIGNIFICANT="--include-non-significant"
VERBOSE="-v"
MAX_TASKS_PER_CHILD="50"

usage() { sed -n '3,38p' "$0"; exit 1; }
log() { printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*"; }
require_cmd() { if ! command -v "$1" >/dev/null 2>&1; then log "ERROR: $1 not found"; exit 3; fi; }
require_dir() { local d="$1"; local name="${2:-$1}"; [[ -d "$d" ]] || { log "ERROR: Required directory not found: $name ($d)"; exit 2; }; }
require_file() { local f="$1"; local name="${2:-$1}"; [[ -f "$f" ]] || { log "ERROR: Required file not found: $name ($f)"; exit 2; }; }

PY_RUNNER="python3"
MAMBA_ENV="${MAMBA_ENV:-caas_global_cancer}"
if command -v micromamba >/dev/null 2>&1; then
  PY_RUNNER="micromamba run -n ${MAMBA_ENV} python"
fi

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--alignment-dir) ALIGNMENT_DIR_REAL="$2"; shift 2;;
    --alignments-toy) ALIGNMENT_DIR_TOY="$2"; shift 2;;
    -t|--tree) TREE_FILE="$2"; shift 2;;
    -o|--output-dir) OUTPUT_DIR="$2"; shift 2;;
    --toy) TOY="$2"; shift 2;;
    --tag) TAG="$2"; shift 2;;
    --taxid-mapping) TAXID_MAPPING="$2"; shift 2;;
    --caas-metadata) CAAS_METADATA="$2"; shift 2;;
    --trait-file) TRAIT_FILE="$2"; shift 2;;
    --ensembl-genes-file) ENSEMBL_GENES_FILE="$2"; shift 2;;
    --asr-mode) ASR_MODE="$2"; shift 2;;
    --asr-model) ASR_MODEL="$2"; shift 2;;
    --asr-cache-dir) ASR_CACHE_DIR="$2"; shift 2;;
    --convergence-mode) CONVERGENCE_MODE="$2"; shift 2;;
    --posterior-threshold) POSTERIOR_THRESHOLD="$2"; shift 2;;
    --run-diagnostics) RUN_DIAGNOSTICS="--run-diagnostics"; shift;;
    --skip-gene-lists) SKIP_GENE_LISTS="--skip-gene-lists"; shift;;
    --include-non-significant) INCLUDE_NON_SIGNIFICANT="--include-non-significant"; shift;;
    --threads) THREADS="$2"; shift 2;;
    --workers) WORKERS="$2"; shift 2;;
    --max-tasks-per-child) MAX_TASKS_PER_CHILD="$2"; shift 2;;
    -v | --verbose) VERBOSE="--verbose"; shift;;
    -h | --help) usage;;
    *) log "Unknown option: $1"; usage;;
  esac
done

require_cmd ${PY_RUNNER%% *}

# Handle TOY mode and alignment selection
toy_lc="$(tr '[:upper:]' '[:lower:]' <<< "${TOY}")"
if [[ "${toy_lc}" == "true" ]]; then
  ALIGNMENT_DIR="${ALIGNMENT_DIR_TOY}"
  if [[ -z "${ALIGNMENT_DIR}" || ! -d "${ALIGNMENT_DIR}" ]]; then
    log "WARNING: --toy true but toy alignments not found; falling back to real alignments: ${ALIGNMENT_DIR_REAL}"
    ALIGNMENT_DIR="${ALIGNMENT_DIR_REAL}"
  fi
  [[ -z "$TAG" ]] && TAG="_toy"
else
  ALIGNMENT_DIR="${ALIGNMENT_DIR_REAL}"
fi

# Validate inputs
require_dir "$ALIGNMENT_DIR" "alignment-dir"
require_file "$TREE_FILE" "tree"
require_file "$CAAS_METADATA" "caas-metadata"
require_file "$TRAIT_FILE" "trait-file"
[[ -n "$ASR_CACHE_DIR" ]] && mkdir -p "$ASR_CACHE_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$SCRIPT_DIR/results/${TRAIT}${TAG}_${STAMP}"
mkdir -p "$OUTPUT_DIR"

log "Starting CAAS disambiguation pipeline"
log "TOY mode: ${TOY}  Tag: ${TAG}"
log "Alignment dir (selected): $ALIGNMENT_DIR"
log "Alignment dir (real): $ALIGNMENT_DIR_REAL"
log "Alignment dir (toy): $ALIGNMENT_DIR_TOY"
log "Tree: $TREE_FILE"
log "CAAS metadata: $CAAS_METADATA"
log "Trait file: $TRAIT_FILE"
log "Ensembl genes file: $ENSEMBL_GENES_FILE"
log "Output dir: $OUTPUT_DIR"
log "ASR mode: $ASR_MODE"
log "Convergence mode: $CONVERGENCE_MODE"
log "Workers: ${WORKERS:-auto}"

cmd="${PY_RUNNER} ${SCRIPT_DIR}/disambiguation_main.py"
cmd="$cmd --alignment-dir \"$ALIGNMENT_DIR\""
cmd="$cmd --tree \"$TREE_FILE\""
cmd="$cmd --caas-metadata \"$CAAS_METADATA\""
cmd="$cmd --trait-file \"$TRAIT_FILE\""
cmd="$cmd --ensembl-genes-file \"$ENSEMBL_GENES_FILE\""
cmd="$cmd --output-dir \"$OUTPUT_DIR\""
cmd="$cmd --asr-mode $ASR_MODE"
cmd="$cmd --asr-model \"$ASR_MODEL\""
cmd="$cmd --convergence-mode \"$CONVERGENCE_MODE\""
cmd="$cmd --posterior-threshold $POSTERIOR_THRESHOLD"
cmd="$cmd --threads $THREADS"

[[ -n "$WORKERS" ]] && cmd="$cmd --workers $WORKERS"
[[ -n "$TAXID_MAPPING" ]] && cmd="$cmd --taxid-mapping \"$TAXID_MAPPING\""
[[ -n "$ASR_CACHE_DIR" ]] && cmd="$cmd --asr-cache-dir \"$ASR_CACHE_DIR\""
[[ -n "$RUN_DIAGNOSTICS" ]] && cmd="$cmd $RUN_DIAGNOSTICS"
[[ -n "$SKIP_GENE_LISTS" ]] && cmd="$cmd $SKIP_GENE_LISTS"
[[ -n "$RUN_DIAGNOSTICS" ]] && cmd="$cmd $RUN_DIAGNOSTICS"
[[ -n "$INCLUDE_NON_SIGNIFICANT" ]] && cmd="$cmd $INCLUDE_NON_SIGNIFICANT"
[[ -n "$VERBOSE" ]] && cmd="$cmd $VERBOSE"
[[ -n "$MAX_TASKS_PER_CHILD" ]] && cmd="$cmd --max-tasks-per-child $MAX_TASKS_PER_CHILD"

log "Command: $cmd"
eval "$cmd"

log "Aggregation pipeline completed"
log "Results in: $OUTPUT_DIR"
