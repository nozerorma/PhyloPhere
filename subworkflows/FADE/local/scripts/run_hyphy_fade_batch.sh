#!/usr/bin/env bash
# run_hyphy_fade_batch.sh
# ───────────────────────────────────────────────────────────────────────────
# Run HyPhy FADE on a batch of genes using bash job control.
# Each gene in the manifest is launched as a background subshell; up to
# --workers run concurrently. Individual gene failures are TOLERATED (logged
# and skipped) so the overall batch task always succeeds.
#
# Manifest format (tab-separated, one gene per line):
#   gene_id <TAB> fasta_filename <TAB> annotated_tree_filename
#
# Files are staged by Nextflow into:
#   fastas/<fasta_filename>
#   trees/<annotated_tree_filename>
# ───────────────────────────────────────────────────────────────────────────
set -uo pipefail  # NOT -e: individual gene failures must not abort the batch

batch_id=""
manifest=""
direction=""
workers="1"
runner_mode=""
model="LG"
model_file_arg=""
method="Variational-Bayes"
grid="20"
concentration="0.5"
mcmc_chains="5"
mcmc_chain_length="2000000"
mcmc_burn_in="1000000"
mcmc_samples="1000"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)           batch_id="$2";          shift 2 ;;
        --manifest)           manifest="$2";           shift 2 ;;
        --direction)          direction="$2";          shift 2 ;;
        --workers)            workers="$2";            shift 2 ;;
        --runner-mode)        runner_mode="$2";        shift 2 ;;
        --model)              model="$2";              shift 2 ;;
        --model-file-arg)     model_file_arg="$2";     shift 2 ;;
        --method)             method="$2";             shift 2 ;;
        --grid)               grid="$2";               shift 2 ;;
        --concentration)      concentration="$2";      shift 2 ;;
        --mcmc-chains)        mcmc_chains="$2";        shift 2 ;;
        --mcmc-chain-length)  mcmc_chain_length="$2";  shift 2 ;;
        --mcmc-burn-in)       mcmc_burn_in="$2";       shift 2 ;;
        --mcmc-samples)       mcmc_samples="$2";       shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$direction" || -z "$runner_mode" ]]; then
    echo "Missing required arguments for FADE batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value: $workers" >&2
    exit 1
fi

# ── Build base HyPhy command ─────────────────────────────────────────────────
declare -a base_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("/usr/local/bin/_entrypoint.sh" "hyphy" "fade")
else
    base_cmd=("hyphy" "fade")
fi

# MCMC args only needed when not using Variational-Bayes
declare -a mcmc_args=()
if [[ "$method" != "Variational-Bayes" ]]; then
    mcmc_args=(
        "--chains"       "$mcmc_chains"
        "--chain-length" "$mcmc_chain_length"
        "--burn-in"      "$mcmc_burn_in"
        "--samples"      "$mcmc_samples"
    )
fi

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched FADE task $batch_id (direction: $direction)"
echo "Genes in batch: $gene_count"
echo "Concurrent workers: $workers"

# ── Job-control helpers ──────────────────────────────────────────────────────
# Note: unlike the CT batch scripts, failures are non-fatal here.

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true   # swallow individual gene failures
    done
}

wait_for_all() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
        wait -n 2>/dev/null || true
    done
}

# ── Process each gene ────────────────────────────────────────────────────────
idx=0
while IFS=$'\t' read -r gene_id fasta_name tree_name; do
    [[ -z "${gene_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[FADE_BATCHED] Launching ${gene_id} (${direction}) ($idx/$gene_count)"

    output_json="${gene_id}.${direction}.FADE.json"
    fasta_path="fastas/${fasta_name}"
    tree_path="trees/${tree_name}"

    (
        # model_file_arg may be empty or "--model-file lg.dat"; pass as separate
        # tokens only when non-empty to avoid quoting issues.
        declare -a cmd=(
            "${base_cmd[@]}"
            "--alignment"               "$fasta_path"
            "--tree"                    "$tree_path"
            "--branches"                "Foreground"
            "--model"                   "$model"
            "--method"                  "$method"
            "--grid"                    "$grid"
            "--concentration_parameter" "$concentration"
            "--output"                  "$output_json"
        )
        if [[ -n "$model_file_arg" ]]; then
            # model_file_arg is e.g. "--model-file lg.dat" — split into two tokens
            read -r -a _mfa <<< "$model_file_arg"
            cmd+=("${_mfa[@]}")
        fi
        cmd+=("${mcmc_args[@]}")

        "${cmd[@]}" \
            || echo "[FADE_BATCHED] FADE failed for ${gene_id} (${direction}), skipping"

        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "$output_json" ] || rm -f "$output_json"
        echo "[FADE_BATCHED] Completed ${gene_id} (${direction})"
    ) &

done < "$manifest"

wait_for_all
echo "[FADE_BATCHED] Batch $batch_id finished."
