#!/usr/bin/env bash
set -euo pipefail

batch_id=""
manifest=""
caas_config=""
workers="1"
ali_format=""
runner_mode=""
ct_bin=""
extra_args_file=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)
            batch_id="$2"
            shift 2
            ;;
        --manifest)
            manifest="$2"
            shift 2
            ;;
        --caas-config)
            caas_config="$2"
            shift 2
            ;;
        --workers)
            workers="$2"
            shift 2
            ;;
        --ali-format)
            ali_format="$2"
            shift 2
            ;;
        --runner-mode)
            runner_mode="$2"
            shift 2
            ;;
        --ct-bin)
            ct_bin="$2"
            shift 2
            ;;
        --extra-args-file)
            extra_args_file="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$caas_config" || -z "$ali_format" || -z "$runner_mode" || -z "$ct_bin" ]]; then
    echo "Missing required arguments for discovery batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value for discovery batch runner: $workers" >&2
    exit 1
fi

declare -a extra_args=()
if [[ -n "$extra_args_file" && -s "$extra_args_file" ]]; then
    # The args file may span multiple lines; tr collapses newlines before word-splitting.
    read -r -a extra_args < <(tr '\n' ' ' < "$extra_args_file"; echo)
fi

declare -a base_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("$ct_bin" "ct" "discovery")
else
    base_cmd=("$ct_bin" "discovery")
fi

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched discovery task $batch_id"
echo "Genes in batch: $gene_count"
echo "Concurrent workers: $workers"

terminate_children() {
    local pids
    pids="$(jobs -pr || true)"
    if [[ -n "$pids" ]]; then
        kill $pids 2>/dev/null || true
        wait $pids 2>/dev/null || true
    fi
}

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        if ! wait -n; then
            echo "[DISCOVERY_BATCHED] A child discovery job failed; stopping batch $batch_id" >&2
            terminate_children
            exit 1
        fi
    done
}

wait_for_all() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
        if ! wait -n; then
            echo "[DISCOVERY_BATCHED] A child discovery job failed; stopping batch $batch_id" >&2
            terminate_children
            exit 1
        fi
    done
}

idx=0
while IFS=$'\t' read -r alignment_id alignment_name; do
    [[ -z "${alignment_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[DISCOVERY_BATCHED] Launching $alignment_id ($idx/$gene_count)"

    alignment_path="alignments/$alignment_name"
    (
        "${base_cmd[@]}" \
            -a "$alignment_path" \
            -t "$caas_config" \
            -o "${alignment_id}.output" \
            --background_output "${alignment_id}.background.tsv" \
            --fmt "$ali_format" \
            "${extra_args[@]}"
        echo "[DISCOVERY_BATCHED] Completed $alignment_id"
    ) &
done <"$manifest"

wait_for_all
