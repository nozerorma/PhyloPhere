#!/usr/bin/env bash
# run_hyphy_molerate_batch.sh
# ───────────────────────────────────────────────────────────────────────────
# Run HyPhy MoleRate on a batch of genes using bash job control.
# Each gene in the manifest is launched as a background subshell; up to
# --workers run concurrently. Individual gene failures are TOLERATED (logged
# and skipped) so the overall batch task always succeeds.
#
# Manifest format (tab-separated, one gene per line):
#   gene_id <TAB> fasta_filename
#
# All genes in a batch share the same --tree and --fg-list (one per direction).
# Files are staged by Nextflow into:
#   fastas/<fasta_filename>
#   tree file at the path given by --tree
#   fg_list file at the path given by --fg-list
# ───────────────────────────────────────────────────────────────────────────
set -uo pipefail  # NOT -e: individual gene failures must not abort the batch

batch_id=""
manifest=""
direction=""
workers="1"
runner_mode=""
tree_path=""
fg_list_path=""
hyphy_analyses="/hyphy-analyses"
model="LG"
model_file_arg=""
rv="None"
rate_classes="4"
labeling_strategy="all-descendants"
branch_tests="Yes"
full_model="Yes"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)            batch_id="$2";           shift 2 ;;
        --manifest)            manifest="$2";            shift 2 ;;
        --direction)           direction="$2";           shift 2 ;;
        --workers)             workers="$2";             shift 2 ;;
        --runner-mode)         runner_mode="$2";         shift 2 ;;
        --tree)                tree_path="$2";           shift 2 ;;
        --fg-list)             fg_list_path="$2";        shift 2 ;;
        --hyphy-analyses)      hyphy_analyses="$2";      shift 2 ;;
        --model)               model="$2";               shift 2 ;;
        --model-file-arg)      model_file_arg="$2";      shift 2 ;;
        --rv)                  rv="$2";                  shift 2 ;;
        --rate-classes)        rate_classes="$2";        shift 2 ;;
        --labeling-strategy)   labeling_strategy="$2";  shift 2 ;;
        --branch-tests)        branch_tests="$2";        shift 2 ;;
        --full-model)          full_model="$2";          shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$direction" || -z "$runner_mode" || -z "$tree_path" || -z "$fg_list_path" ]]; then
    echo "Missing required arguments for MoleRate batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value: $workers" >&2
    exit 1
fi

# ── Resolve molerate.bf path ─────────────────────────────────────────────────
MOLERATE_BF="${hyphy_analyses}/molerate/molerate.bf"
if [[ "$runner_mode" != "container" && ! -f "$MOLERATE_BF" ]]; then
    MOLERATE_BF="$(dirname "$(which hyphy)")/../share/hyphy/TemplateBatchFiles/molerate.bf"
fi

# ── Build shared --branches args array from fg_list ──────────────────────────
# Using an array is safe: no word-splitting or injection risk.
declare -a BRANCHES_ARGS=()
while IFS= read -r branch || [[ -n "$branch" ]]; do
    branch="${branch%%$'\r'}"   # strip Windows CR if present
    [[ -z "$branch" ]] && continue
    BRANCHES_ARGS+=("--branches" "$branch")
done < "$fg_list_path"

# ── Build base HyPhy command ─────────────────────────────────────────────────
declare -a base_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("/usr/local/bin/_entrypoint.sh" "hyphy" "$MOLERATE_BF")
else
    base_cmd=("hyphy" "$MOLERATE_BF")
fi

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched MoleRate task $batch_id (direction: $direction)"
echo "Genes in batch: $gene_count"
echo "Concurrent workers: $workers"

# ── Job-control helpers ──────────────────────────────────────────────────────
# Individual gene failures are non-fatal.

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

wait_for_all() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
        wait -n 2>/dev/null || true
    done
}

# ── Process each gene ────────────────────────────────────────────────────────
idx=0
while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "${gene_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[MOLERATE_BATCHED] Launching ${gene_id} (${direction}) ($idx/$gene_count)"

    output_json="${gene_id}.${direction}.molerate.json"
    fasta_path="fastas/${fasta_name}"

    # Capture the arrays for the subshell (bash subshells inherit variables but
    # not indexed arrays declared after the subshell forks; capture them here).
    declare -a _branches_snapshot=("${BRANCHES_ARGS[@]}")
    declare -a _base_snapshot=("${base_cmd[@]}")

    (
        declare -a cmd=(
            "${_base_snapshot[@]}"
            "--alignment"              "$fasta_path"
            "--tree"                   "$tree_path"
            "--model"                  "$model"
            "--rv"                     "$rv"
            "--rate-classes"           "$rate_classes"
            "--labeling-strategy"      "$labeling_strategy"
            "--branch-level-analysis"  "$branch_tests"
            "--full-model"             "$full_model"
            "--output"                 "$output_json"
        )
        if [[ -n "$model_file_arg" ]]; then
            # model_file_arg is e.g. "--model-file lg.dat" — split into two tokens
            read -r -a _mfa <<< "$model_file_arg"
            cmd+=("${_mfa[@]}")
        fi
        cmd+=("${_branches_snapshot[@]}")

        "${cmd[@]}" \
            || echo "[MOLERATE_BATCHED] MoleRate failed for ${gene_id} (${direction}), skipping"

        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "$output_json" ] || rm -f "$output_json"
        echo "[MOLERATE_BATCHED] Completed ${gene_id} (${direction})"
    ) &

done < "$manifest"

wait_for_all
echo "[MOLERATE_BATCHED] Batch $batch_id finished."
