#!/bin/bash
# Reproduce the Tier 0 gate end to end.
#
#   validation/tier0/reproduce.sh [OUT_DIR] [N_REPLICATES] [N_GENES]
#
# Defaults reproduce the reference run in validation/tier0/gate_result/
# (10 replicates/cell, 120 genes). Needs the primate tree fixture
# (validation/fixtures/tier0/README.md) and the phylophere conda env.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

OUT="${1:-validation/runs/tier0_gate}"
NREP="${2:-10}"
NGENES="${3:-120}"
LAMBDAS="${LAMBDAS:-0,0.5,1}"
JOBS="${JOBS:-2}"

echo "=== stage + run: $OUT ($NREP reps/cell, $NGENES genes, lambda=$LAMBDAS, jobs=$JOBS) ==="
python -m validation.tier0.run_replicates --out "$OUT" \
    --archetypes binary,rate,continuous --sets null,power --trees primate_half \
    --lambdas "$LAMBDAS" --n-replicates "$NREP" --n-genes "$NGENES" \
    --n-sites 350 --n-pairs 4 --run --jobs "$JOBS"

echo "=== score ==="
python -m validation.harness.cli score --run "$OUT" --json "$OUT/score.json"

echo "=== figures ==="
python -m validation.tier0.plots --run "$OUT" --out "$OUT/figures"

echo
echo "score.json : $OUT/score.json"
echo "figures    : $OUT/figures/"
