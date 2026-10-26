#!/bin/bash
# Profile POSITION_REGISTRY (posenrich_registry.py) on the staged inputs of one Nextflow task.
# Run inside an allocation (salloc / srun --pty), never on the login node.
# Usage: profile_posenrich_registry.sh TASK_WORKDIR PROFILE_OUT
#   TASK_WORKDIR  the POSITION_REGISTRY work dir (holds the staged inputs as symlinks)
#   PROFILE_OUT   new directory below $HOME/scratch/0.Phylophere/.tmp (created if missing)
# Environment: PYTHON selects the interpreter that has pandas (default: python3).
set -euo pipefail

# Refuse to run outside a Slurm allocation, as the cluster policy requires
[ -n "${SLURM_JOB_ID:-}" ] || { echo 'Run inside salloc or srun --pty, not on the login node' >&2; exit 2; }

task=$(realpath "$1")
out=$(realpath -m "$2")
tmp_root="$HOME/scratch/0.Phylophere/.tmp"
case "$out" in "$tmp_root"/*) ;; *) echo "PROFILE_OUT must be below $tmp_root" >&2; exit 2;; esac
[ ! -e "$out" ] || { echo "PROFILE_OUT already exists: $out" >&2; exit 2; }
mkdir -p "$out/registry"

# Code under test: the checkout this script belongs to
code=$(realpath "$(dirname "$0")/../..")
registry_py="$code/subworkflows/ENRICHMENT/local/src/posenrich_registry.py"
PYTHON=${PYTHON:-python3}

# Same arguments as the task's .command.sh, resolved against the task work dir
cd "$task"
args=(
    --position_scores position_scores.tsv
    --postproc filtered_discovery.tsv
    --null_universe null_universe_summary.tsv
    --map_dir MAP
    --cleaned_background cleaned_background_main.txt
    --gene_ensembl_file ensembl_genes.output
    --domain_variability_file domain_variability.tsv.gz
    --ucr_positions_file ucr_positions.tsv.gz
    --fubar_sites_file fubar_sites.tsv.gz
    --fade_sites_top_file fade_sites_top.csv
    --fade_sites_bottom_file fade_sites_bottom.csv
    --custom_marker_file pp2a_custom.tsv
    --output_dir "$out/registry"
    --gmt_dir "$out/registry/gmts"
)

# Wall time and peak memory of the profiled run. cProfile adds overhead, so the
# absolute times are inflated; use the profile for relative weights only.
export TMPDIR="$out"
status=0
if [ -x /usr/bin/time ]; then
    /usr/bin/time -v -o "$out/time.txt" "$PYTHON" -m cProfile -o "$out/registry.prof" "$registry_py" "${args[@]}" \
        > "$out/run.log" 2>&1 || status=$?
else
    start=$(date +%s)
    "$PYTHON" -m cProfile -o "$out/registry.prof" "$registry_py" "${args[@]}" > "$out/run.log" 2>&1 || status=$?
    echo "wall_s=$(( $(date +%s) - start ))" > "$out/time.txt"
fi
# A failed run still leaves run.log for diagnosis; report it instead of hiding it
[ "$status" -eq 0 ] || { echo "registry exited with status $status; see $out/run.log" >&2; exit "$status"; }

# Human-readable summary: top functions by cumulative and by own time
"$PYTHON" -c "
import pstats, sys
prof = sys.argv[1]
with open(sys.argv[2], 'w') as f:
    for key in ('cumulative', 'tottime'):
        stats = pstats.Stats(prof, stream=f)
        stats.sort_stats(key).print_stats(40)
" "$out/registry.prof" "$out/profile_summary.txt"

echo "Profile written to $out (see profile_summary.txt, time.txt, run.log)"
