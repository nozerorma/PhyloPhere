#!/bin/bash
# Run inside an allocation/job, with the synchronized checkout on PATH/configured.
# Usage: prepare_snapshot.sh CHECKOUT SCRATCH_OUT [Nextflow arguments for wp6_validation.nf: --wp6_stats --wp6_universe --wp6_perms]
set -euo pipefail
checkout=$(realpath "$1")
output=$(realpath -m "$2")
shift 2
case "$output" in "$HOME"/scratch/*) ;; *) echo 'Output must be below $HOME/scratch' >&2; exit 2;; esac
[ ! -e "$output/code" ] || { echo "Snapshot already exists: $output/code" >&2; exit 2; }
mkdir -p "$output/code" "$output/tmp"
export TMPDIR="$output/tmp"
git -C "$checkout" archive HEAD | tar -x -C "$output/code"
cat > "$output/code/wp6_validation.nf" <<'NF'
nextflow.enable.dsl=2
include { FCS_RER } from './subworkflows/ENRICHMENT/fcs_inputs.nf'
workflow {
 FCS_RER(Channel.value(file(params.wp6_stats)),Channel.value(file(params.wp6_universe)),Channel.value(file(params.wp6_perms)))
 FCS_RER.out.enrich_file.subscribe { f -> println "WP6_RESULT=${f}" }
}
NF
cd "$output/code"
git -C "$checkout" rev-parse HEAD > "$output/commit.txt"
nextflow -log "$output/nextflow.log" run wp6_validation.nf -work-dir "$output/work" \
 -with-trace "$output/trace.tsv" -with-report "$output/execution.html" -with-timeline "$output/timeline.html" "$@"
