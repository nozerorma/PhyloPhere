#!/usr/bin/env bash
set -euo pipefail

batch_id=""
manifest=""
caas_config=""
resampled_path=""
workers="1"
ali_format=""
ct_bin=""
export_groups="0"
export_perm_discovery="0"
export_b0="0"
export_universe="0"
export_universe_strict="0"
extra_args_file=""
stall_timeout="1800"

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
        --ct-bin)
            ct_bin="$2"
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
        --export-b0)
            export_b0="$2"
            shift 2
            ;;
        --export-universe)
            export_universe="$2"
            shift 2
            ;;
        --export-universe-strict)
            export_universe_strict="$2"
            shift 2
            ;;
        --extra-args-file)
            extra_args_file="$2"
            shift 2
            ;;
        --stall-timeout)
            stall_timeout="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$caas_config" || -z "$resampled_path" || -z "$ali_format" || -z "$ct_bin" ]]; then
    echo "Missing required arguments for perm-replay batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value for perm-replay batch runner: $workers" >&2
    exit 1
fi

declare -a extra_args=()
if [[ -n "$extra_args_file" && -s "$extra_args_file" ]]; then
    # The args file may span multiple lines; tr collapses newlines before word-splitting.
    read -r -a extra_args < <(tr '\n' ' ' < "$extra_args_file"; echo)
fi

# ct_bin may be a single path (local mode) or a space-separated command
# prefix such as "/usr/local/bin/_entrypoint.sh <path-to-ct>" (container mode).
declare -a ct_bin_arr
prefix_file=".ct_prefix_${batch_id}.argv"
python3 - "$ct_bin" > "$prefix_file" <<'PHYLOPHERE_PREFIX'
import shlex, sys
parts = shlex.split(sys.argv[1])
if not parts:
    raise SystemExit('Empty CT command prefix')
sys.stdout.buffer.write(b'\0'.join(p.encode() for p in parts) + b'\0')
PHYLOPHERE_PREFIX
mapfile -d '' -t ct_bin_arr < "$prefix_file"
rm -f "$prefix_file"
declare -a base_cmd=("${ct_bin_arr[@]}" "perm-replay")

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched perm-replay task $batch_id"
echo "Genes in batch: $gene_count"
echo "Concurrent workers: $workers (each parses the labelings once)"

terminate_children() {
    local pids
    pids="$(jobs -pr || true)"
    if [[ -n "$pids" ]]; then
        kill $pids 2>/dev/null || true
        wait $pids 2>/dev/null || true
    fi
}

# Kills $1 if it burns zero CPU time for $stall_timeout seconds straight.
# Guards against workers wedged in a filesystem/import stall (0% CPU,
# indefinite openat() block on the shared conda env) rather than genuinely
# slow computation, which keeps accruing CPU time and is left alone.
watchdog_guard() {
    local pid="$1" poll=30 last_cpu=-1 same=0
    local max_same=$(( stall_timeout / poll ))
    while kill -0 "$pid" 2>/dev/null; do
        sleep "$poll"
        local cpu
        cpu="$(awk '{print $14+$15}' "/proc/$pid/stat" 2>/dev/null)" || break
        [[ -z "$cpu" ]] && break
        if [[ "$cpu" == "$last_cpu" ]]; then
            same=$((same + 1))
            if [[ "$same" -ge "$max_same" ]]; then
                echo "[perm-replay batch] Worker pid $pid stalled at 0% CPU for ${stall_timeout}s; killing" >&2
                kill -TERM "$pid" 2>/dev/null || true
                sleep 5
                kill -KILL "$pid" 2>/dev/null || true
                break
            fi
        else
            same=0
        fi
        last_cpu="$cpu"
    done
}

wait_for_all() {
    local pid status
    # Completed jobs still need waiting: jobs -pr omits failures that finished
    # before collection. Keep every launched PID, including fast failures.
    for pid in "${worker_pids[@]}"; do
        status=0
        wait "$pid" || status=$?
        if [[ "$status" -ne 0 ]]; then
            echo "[perm-replay batch] A child failed (exit $status); stopping $batch_id" >&2
            terminate_children
            exit "$status"
        fi
    done
}

# Each of the `workers` processes parses the labelings once and then takes genes from the manifest (a gene is taken
# by creating its directory under claims/, so the genes spread over the processes as they finish). The export paths are
# templates: ct replaces {id} with the gene id.
claim_dir="claims_${batch_id}"
rm -rf "$claim_dir"
mkdir -p "$claim_dir"

declare -a cmd=(
    "${base_cmd[@]}"
    --manifest "$manifest"
    --alignment_dir alignments
    --claim_dir "$claim_dir"
    -t "$caas_config"
    -s "$resampled_path"
    --fmt "$ali_format"
)

if [[ "$export_groups" == "1" ]]; then
    cmd+=(--export_groups '{id}.perm_replay.groups.output')
fi
if [[ "$export_perm_discovery" == "1" ]]; then
    cmd+=(--export_perm_discovery '{id}.perm_replay.discovery.output')
fi
if [[ "$export_b0" == "1" ]]; then
    # the discovery.tab rows and the tested positions of the b_0 labeling(s); no discovery file when b_0 has no hit
    cmd+=(--export_b0_discovery '{id}.b0.discovery.tsv' --export_b0_background '{id}.b0.background')
fi

if [[ "$export_universe" == "1" ]]; then
    # the positional universe of the b_0 hypotheses (positions tested by b_0 and eligible under some pair orientation)
    cmd+=(--export_b0_universe '{id}.b0.universe')
fi

if [[ "$export_universe_strict" == "1" ]]; then
    # the same universe with no conserved pair allowed (the strict analysis)
    cmd+=(--export_b0_universe_strict '{id}.b0.universe_strict')
fi

declare -a worker_pids=()
trap terminate_children INT TERM
n_proc=$(( gene_count < workers ? gene_count : workers ))
for (( p = 1; p <= n_proc; p++ )); do
    (
        "${cmd[@]}" "${extra_args[@]}" &
        worker_pid=$!
        ( watchdog_guard "$worker_pid" ) >/dev/null 2>&1 &
        watchdog_pid=$!
        disown "$watchdog_pid" 2>/dev/null || true
        worker_status=0
        wait "$worker_pid" || worker_status=$?
        kill "$watchdog_pid" 2>/dev/null || true
        exit "$worker_status"
    ) </dev/null &
    worker_pids+=("$!")
done

wait_for_all
rm -rf "$claim_dir"
