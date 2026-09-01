#!/bin/bash
# Tier 1 — first real PEPC pipeline run (Task 2).
#
#   validation/tier1/run_pepc.sh 0.75      # min_divergent_fraction sweep value
#   validation/tier1/run_pepc.sh 1.0
#
# Module set (build_project.py): contrast_selection + CAAS/CAAP(disc,resample,
# bootstrap,permulation) + CT_DISAMBIGUATION(+POSTPROC + ASR robustness) + FADE
# + SCORING. RER / accumulation / enrichment / VEP / reporting OFF.
#
# Fixture must be built first:
#   micromamba run -n phylophere-fixtures python3 \
#       validation/fixtures/tier1/pepc/build.py
set -Eeuo pipefail

MDF="${1:?usage: run_pepc.sh <min_divergent_fraction>  e.g. 0.75 or 1.0}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FIX="${REPO_DIR}/validation/fixtures/tier1/pepc"
GEN="${REPO_DIR}/validation/runs/tier1/_generated"
OUT="${REPO_DIR}/validation/runs/tier1/pepc/mdf${MDF}"

for f in "${FIX}/align/PEPC.fasta" "${FIX}/tree.nwk" "${FIX}/my_traits.tsv"; do
    [ -e "$f" ] || { echo "missing fixture file: $f  (run pepc/build.py)" >&2; exit 1; }
done

# (re)render run_phenotype_single.sh + the env preamble for this mdf
micromamba run -n phylophere python3 "${REPO_DIR}/validation/tier1/build_project.py" \
    --mdf "${MDF}" --out "${GEN}"

mkdir -p "${OUT}"
set -a
source "${GEN}/tier0_env.sh"          # render() hardcodes this filename
ALI_DIR="${FIX}/align"
ALI_FORMAT="fasta"
ALI_SP_NAMES="${FIX}/ali_sp_names.txt"
INPUT_TAX_ID="${FIX}/taxid.tsv"
GENE_ENSEMBL_FILE="${FIX}/gene_ensembl.tsv"
TREE_FILE="${FIX}/tree.nwk"
TRAIT_FILE="${FIX}/my_traits.tsv"
SIMPLE_TRAIT_FILE="${FIX}/my_traits.tsv"
MIN_DIVERGENT_FRACTION="${MDF}"
SEED="1998"
CAAS_OUTBASE="${OUT}/out"
WORK_BASE="${OUT}/work"
ASR_CACHE_DIR="${OUT}/asr_cache"
set +a
mkdir -p "$CAAS_OUTBASE" "$WORK_BASE" "$ASR_CACHE_DIR"

echo "=== PEPC Tier 1 run  min_divergent_fraction=${MDF} ==="
echo "  out : ${OUT}"
#            CLASS TRAIT SECONDARY N_TRAIT C_TRAIT PRUNE PRUNE_SEC DISCRETE TRAIT_TYPE
exec bash "${REPO_DIR}/run_phenotype_single.sh" 2 c4 "" "" "" "" "" "" ordinal
