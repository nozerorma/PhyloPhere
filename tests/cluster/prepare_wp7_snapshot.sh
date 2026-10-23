#!/bin/bash
# Run inside an allocation: prepare_wp7_snapshot.sh CHECKOUT SCRATCH_OUT SAVED_RUN [Nextflow arguments]
set -euo pipefail
checkout=$(realpath "$1")
output=$(realpath -m "$2")
saved=$(realpath "$3")
shift 3
case "$output" in "$HOME"/scratch/*) ;; *) echo 'Output must be below $HOME/scratch' >&2; exit 2;; esac
[ ! -e "$output/code" ] || { echo "Snapshot already exists: $output/code" >&2; exit 2; }
mkdir -p "$output/code" "$output/tmp"
export TMPDIR="$output/tmp"
git -C "$checkout" archive HEAD | tar -x -C "$output/code"
sed "s|from '../|from './|g" "$output/code/tests/wp7_validation.nf" > "$output/code/wp7_validation.nf"
git -C "$checkout" rev-parse HEAD > "$output/commit.txt"
cd "$output/code"
nextflow -log "$output/nextflow.log" run wp7_validation.nf -work-dir "$output/work" \
 --wp7_saved "$saved" --outdir "$output/results" -with-trace "$output/trace.tsv" -with-report "$output/execution.html" -with-timeline "$output/timeline.html" "$@"
