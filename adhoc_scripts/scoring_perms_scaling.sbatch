#!/bin/bash
# Compare validation-only (normal) with the whole-memory compatibility export on a canonical store.
# Submit inside the cluster scheduler; SOURCE is a caas_permulation directory, MODE=validate|export.
#SBATCH --job-name=scoring_perms_scaling
#SBATCH --partition=std-cpu
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o /homes/users/mramon/scratch/0.Phylophere/.tmp/scoring_perms_scaling/slurm-%A_%a.out
#SBATCH -e /homes/users/mramon/scratch/0.Phylophere/.tmp/scoring_perms_scaling/slurm-%A_%a.err

set -euo pipefail
ROOT="$HOME/scratch/0.Phylophere"
: "${SOURCE:?Set SOURCE to a canonical caas_permulation directory}"
MODE="${MODE:-validate}"
OUT="$ROOT/.tmp/scoring_perms_scaling/${SLURM_JOB_ID}_${MODE}"
mkdir -p "$OUT/tmp"
export TMPDIR="$OUT/tmp"
flags=(--validate-only)
if [ "$MODE" = export ]; then flags=(--output "$OUT/caas_perms.rds"); fi
if [ "$MODE" != validate ] && [ "$MODE" != export ]; then echo 'MODE must be validate or export' >&2; exit 2; fi
"$HOME/micromamba/envs/phylophere/bin/Rscript" "$ROOT/subworkflows/SCORING/local/src/scoring_caas_perms.R" \
 --blocks "$SOURCE"/blocks/blk_*/matrices.rds --cycles "$SOURCE/cycle_roster.txt" \
 --metadata "$OUT/null_metadata.json" "${flags[@]}"
