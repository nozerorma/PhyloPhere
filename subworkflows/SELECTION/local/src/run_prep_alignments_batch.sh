#!/usr/bin/env bash
# Compatibility entrypoint; Python owns validation, bounded workers and failure propagation.
set -euo pipefail
batch_id="" manifest="" workers="" runner_mode="" tree="" local_dir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id) batch_id="$2"; shift 2 ;;
        --manifest) manifest="$2"; shift 2 ;;
        --workers) workers="$2"; shift 2 ;;
        --runner-mode) runner_mode="$2"; shift 2 ;;
        --tree) tree="$2"; shift 2 ;;
        --local-dir) local_dir="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$batch_id" && -n "$manifest" && -n "$tree" && -n "$local_dir" ]] || { echo "Missing required arguments" >&2; exit 1; }
[[ "$workers" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid workers" >&2; exit 1; }
case "$runner_mode" in
    container) py=(/usr/local/bin/_entrypoint.sh python) ;;
    local) py=(python) ;;
    *) echo "Invalid runner mode" >&2; exit 1 ;;
esac
exec "${py[@]}" -B "$local_dir/prepare_alignments.py" --manifest "$manifest" --workers "$workers" --tree "$tree"
