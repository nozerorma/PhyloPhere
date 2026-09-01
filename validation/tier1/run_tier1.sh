#!/bin/bash
# Tier 1 — drive one pipeline run over a fixture (CAAS + FADE, no RER).
#
#   run_tier1.sh <fixture> <traitcol> <trait_type> <mdf> [asr_model] [tag] [qtail]
#
#     fixture     dir name under validation/fixtures/tier1/
#     traitcol    column of that fixture's my_traits.tsv to use as --traitname
#     trait_type  ordinal | continuous
#     mdf         min_divergent_fraction (0.75 / 1.0)
#     asr_model   lg (default) | jtt | wag       -- ASR-robustness control
#     tag         optional extra suffix for the run dir (e.g. wrongfg)
#     discrete    continuous only: discrete_method (default quintile -> q20/q80;
#                 also quartile / decile / median_sd / parameterized). Setting a
#                 named method overrides the top/bottom_quantile params. Ordinal
#                 traits ignore this (extreme levels are the cut-offs).
#
# Examples:
#   run_tier1.sh pepc        c4       ordinal    0.75
#   run_tier1.sh pepc        c4       ordinal    0.75 jtt
#   run_tier1.sh hb_altitude elev_mid continuous 0.75 lg "" quartile
#   run_tier1.sh hb_altitude altitude ordinal    0.75
#
# Fixture must be built first (validation/fixtures/tier1/<fixture>/build.py).
set -Eeuo pipefail

FIXNAME="${1:?fixture}"; TRAITCOL="${2:?traitcol}"; TTYPE="${3:?trait_type}"
MDF="${4:?min_divergent_fraction}"; ASR="${5:-lg}"; TAG="${6:-}"

if [ "$TTYPE" = continuous ]; then
    DMETHOD="${7:-quintile}"
else
    DMETHOD=""                             # ordinal: unused
fi

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FIX="${REPO_DIR}/validation/fixtures/tier1/${FIXNAME}"
GEN="${REPO_DIR}/validation/runs/tier1/_generated"
SUF="${TRAITCOL}_mdf${MDF}"; [ "$ASR" != lg ] && SUF="${SUF}_${ASR}"; [ -n "$TAG" ] && SUF="${SUF}_${TAG}"
OUT="${REPO_DIR}/validation/runs/tier1/${FIXNAME}/${SUF}"

[ -d "${FIX}/align" ] || { echo "no ${FIX}/align — run ${FIXNAME}/build.py" >&2; exit 1; }
[ -f "${FIX}/my_traits.tsv" ] || { echo "no ${FIX}/my_traits.tsv" >&2; exit 1; }
head -1 "${FIX}/my_traits.tsv" | tr '\t' '\n' | grep -qx "${TRAITCOL}" \
    || { echo "column '${TRAITCOL}' not in ${FIX}/my_traits.tsv" >&2; exit 1; }

micromamba run -n phylophere python3 "${REPO_DIR}/validation/tier1/build_project.py" \
    --trait "${TRAITCOL}" --trait-type "${TTYPE}" --mdf "${MDF}" --asr-model "${ASR}" \
    --clade "tier1_${FIXNAME}" --out "${GEN}"

mkdir -p "${OUT}"
set -a
source "${GEN}/tier0_env.sh"
ALI_DIR="${FIX}/align"
ALI_FORMAT="fasta"
ALI_SP_NAMES="${FIX}/ali_sp_names.txt"
INPUT_TAX_ID="${FIX}/taxid.tsv"
GENE_ENSEMBL_FILE="${FIX}/gene_ensembl.tsv"
TREE_FILE="${FIX}/tree.nwk"
TRAIT_FILE="${FIX}/my_traits.tsv"
SIMPLE_TRAIT_FILE="${FIX}/my_traits.tsv"
MIN_DIVERGENT_FRACTION="${MDF}"
ASR_MODEL="${ASR}"
DISCRETE_METHOD="${DMETHOD}"
SEED="1998"
CAAS_OUTBASE="${OUT}/out"
WORK_BASE="${OUT}/work"
ASR_CACHE_DIR="${OUT}/asr_cache"
set +a
mkdir -p "$CAAS_OUTBASE" "$WORK_BASE" "$ASR_CACHE_DIR"

echo "=== Tier 1  ${FIXNAME} :: ${TRAITCOL} (${TTYPE})  mdf=${MDF}  asr=${ASR}  tag='${TAG}' ==="
echo "  out    : ${OUT}"
[ "$TTYPE" = continuous ] && echo "  discret: ${DMETHOD}"
#            CLASS TRAIT SECONDARY N_TRAIT C_TRAIT PRUNE PRUNE_SEC DISCRETE      TRAIT_TYPE
exec bash "${REPO_DIR}/run_phenotype_single.sh" 2 "${TRAITCOL}" "" "" "" "" "" "${DMETHOD}" "${TTYPE}"
