#!/usr/bin/env bash
set -euo pipefail

batch_id=""
manifest=""
caas_config=""
resampled_path=""
workers="1"
ali_format=""
runner_mode=""
ct_bin=""
progress_log="0"
export_groups="0"
export_perm_discovery="0"
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
        --resampled-path)
            resampled_path="$2"
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
        --progress-log)
            progress_log="$2"
            shift 2
            ;;
        --export-groups)
            export_groups="$2"
            shift 2
            ;;
        --export-perm-discovery)
            export_perm_discovery="$2"
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

if [[ -z "$batch_id" || -z "$manifest" || -z "$caas_config" || -z "$resampled_path" || -z "$ali_format" || -z "$runner_mode" || -z "$ct_bin" ]]; then
    echo "Missing required arguments for bootstrap batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value for bootstrap batch runner: $workers" >&2
    exit 1
fi

declare -a extra_args=()
if [[ -n "$extra_args_file" && -s "$extra_args_file" ]]; then
    # The args file may span multiple lines; tr collapses newlines before word-splitting.
    read -r -a extra_args < <(tr '\n' ' ' < "$extra_args_file"; echo)
fi

declare -a base_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("$ct_bin" "ct" "bootstrap")
else
    base_cmd=("$ct_bin" "bootstrap")
fi

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched bootstrap task $batch_id"
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
            echo "[BOOTSTRAP_BATCHED] A child bootstrap job failed; stopping batch $batch_id" >&2
            terminate_children
            exit 1
        fi
    done
}

wait_for_all() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
        if ! wait -n; then
            echo "[BOOTSTRAP_BATCHED] A child bootstrap job failed; stopping batch $batch_id" >&2
            terminate_children
            exit 1
        fi
    done
}

idx=0
while IFS=$'\t' read -r alignment_id alignment_name discovery_name; do
    [[ -z "${alignment_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[BOOTSTRAP_BATCHED] Launching $alignment_id ($idx/$gene_count)"

    alignment_path="alignments/$alignment_name"

    declare -a cmd=(
        "${base_cmd[@]}"
        -a "$alignment_path"
        -t "$caas_config"
        -s "$resampled_path"
        -o "${alignment_id}.bootstraped.output"
        --fmt "$ali_format"
    )

    if [[ "$discovery_name" != "NO_FILE" ]]; then
        cmd+=(--discovery "discovery/$discovery_name")
    fi
    if [[ "$progress_log" == "1" ]]; then
        cmd+=(--progress_log "${alignment_id}.progress.log")
    fi
    if [[ "$export_groups" == "1" ]]; then
        cmd+=(--export_groups "${alignment_id}.bootstrap.groups.output")
    fi
    if [[ "$export_perm_discovery" == "1" ]]; then
        cmd+=(--export_perm_discovery "${alignment_id}.bootstrap.discovery.output")
    fi

    (
        "${cmd[@]}" "${extra_args[@]}"
        echo "[BOOTSTRAP_BATCHED] Completed $alignment_id"
    ) &
done <"$manifest"

wait_for_all
