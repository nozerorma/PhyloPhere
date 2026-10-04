#!/usr/bin/env bash
# Gate 2.2 review: the toy run of the b_0 core against its baseline, with the three comparison tools.
#
#   review_gate_2_2.sh [--dry-run] <baseline results dir> <new results dir> [report dir]
#
# Runs, in this order, writing one JSON report each into the report dir (default: <new results dir>/gate_2_2):
#   compare_b0.py        <new>                  observed vs the b_0 slice of the same run, checkpoints A-E
#   compare_null.py      --a <baseline> --b <new>   the null tables and the scoring tables, bit for bit
#   compare_contract.py  --a <baseline> --b <new>   the observed contract files
# and exits 1 if any of them fails. It reads large tables: on a cluster run it through srun or sbatch, never on the
# login node (it refuses there unless ALLOW_LOGIN_NODE is set), e.g.
#
#   srun -p high-cpu -c 2 --mem=8G -t 20 bash validation/unification/review_gate_2_2.sh <baseline> <new>
#
# PYTHON selects the interpreter (default python3: the pipeline's environment must provide pandas and numpy).
set -euo pipefail

dry=0
if [[ "${1:-}" == "--dry-run" ]]; then dry=1; shift; fi
if [[ $# -lt 2 || $# -gt 3 ]]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 2
fi
base="$1"; new="$2"; out="${3:-$new/gate_2_2}"
host="${REVIEW_HOSTNAME:-$(hostname)}"
if [[ -z "${SLURM_JOB_ID:-}" && -z "${ALLOW_LOGIN_NODE:-}" && "$host" == correfoc* ]]; then
    echo "Refusing to compare on the login node $host: use srun or sbatch (see the header of this script)." >&2
    exit 2
fi

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${PYTHON:-python3}"
cmds=(
    "$py $here/compare_b0.py --run $new --out $out/compare_b0.json"
    "$py $here/compare_null.py --a $base --b $new --extra scoring/position_scores.tsv --extra scoring/gene_scores.tsv"
    "$py $here/compare_contract.py --a $base --b $new --report $out/compare_contract.json"
)
if [[ $dry -eq 1 ]]; then printf '%s\n' "${cmds[@]}"; exit 0; fi

mkdir -p "$out"
status=0
for cmd in "${cmds[@]}"; do
    echo "=== $cmd"
    name="$(echo "$cmd" | sed -E 's/.*(compare_[a-z0-9]+)\.py.*/\1/')"
    if ! $cmd 2>&1 | tee "$out/$name.log"; then status=1; fi
done
echo "=== reports and logs in $out"
exit $status
