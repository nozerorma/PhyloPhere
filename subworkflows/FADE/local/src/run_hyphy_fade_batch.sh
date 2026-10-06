#!/usr/bin/env bash
# run_hyphy_fade_batch.sh — Run HyPhy FADE on a batch of genes with bash job control.
# PhyloPhere | subworkflows/FADE/local/src/
# =============================================================================
# Called by:  FADE_BATCHED Nextflow process (fade_run.nf → bash run_hyphy_fade_batch.sh ...)
#
# Each gene of the manifest runs as a background subshell, up to --workers at a time.
# A failed gene is logged and skipped, so the batch task itself succeeds; an empty
# JSON is removed so that it is not emitted as an output.
#
# Manifest (tab-separated, one gene per line, written by the process):
#   gene_id <TAB> fasta_filename <TAB> annotated_tree_filename
# The files are staged by Nextflow as fastas/<fasta_filename> and trees/<annotated_tree_filename>.
#
# Args (named flags):
#   --batch-id, --manifest, --direction (top | bottom), --runner-mode (container | local)
#       required; container mode runs hyphy through /usr/local/bin/_entrypoint.sh
#   --workers            concurrent genes (positive integer, default 1)
#   --cpu-per-worker     CPUs of each HyPhy call (default 1)
#   --model, --model-file-arg, --method, --grid, --concentration
#       FADE options (--model-file-arg is split into separate tokens)
#   --mcmc-chains, --mcmc-chain-length, --mcmc-burn-in, --mcmc-samples
#       MCMC options, used only when --method is not Variational-Bayes
# Output: <gene_id>.<direction>.FADE.json in the working directory.
# =============================================================================
set -uo pipefail  # no -e: a failed gene must not abort the batch

batch_id=""
manifest=""
direction=""
workers="1"
cpu_per_worker="1"
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
        --cpu-per-worker)     cpu_per_worker="$2";     shift 2 ;;
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

# ── Base HyPhy command ───────────────────────────────────────────────────────
declare -a base_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("/usr/local/bin/_entrypoint.sh" "hyphy" "fade")
else
    base_cmd=("hyphy" "fade")
fi

# MCMC options are passed only for the MCMC method
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
echo "Concurrent workers: $workers (${cpu_per_worker} CPU(s) per HyPhy call)"

# ── Job-control helpers ──────────────────────────────────────────────────────
# A failed gene is not fatal; the loops swallow the exit status of finished jobs.

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
        # model_file_arg is empty or "--model-file lg.dat"; it is added as separate
        # tokens only when non-empty, to avoid quoting problems.
        declare -a cmd=(
            "${base_cmd[@]}"
            "--alignment"               "$fasta_path"
            "--tree"                    "$tree_path"
            "--branches"                "Foreground"
            "--model"                   "$model"
            "--method"                  "$method"
            "--grid"                    "$grid"
            "--concentration_parameter" "$concentration"
            "--cpu"                     "$cpu_per_worker"
            "--output"                  "$output_json"
        )
        if [[ -n "$model_file_arg" ]]; then
            # split "--model-file lg.dat" into two tokens
            read -r -a _mfa <<< "$model_file_arg"
            cmd+=("${_mfa[@]}")
        fi
        cmd+=("${mcmc_args[@]}")

        # Cap BLAS and OpenMP threads at the per-worker CPUs: --cpu only sets the
        # HyPhy scheduler, while OpenBLAS, OpenMP and MKL read the CPU count of the
        # node unless these variables are set.
        export OMP_NUM_THREADS="${cpu_per_worker}"
        export MKL_NUM_THREADS="${cpu_per_worker}"
        export OPENBLAS_NUM_THREADS="${cpu_per_worker}"
        export BLAS_NUM_THREADS="${cpu_per_worker}"

        "${cmd[@]}" \
            || echo "[FADE_BATCHED] FADE failed for ${gene_id} (${direction}), skipping"

        # Remove 0-byte JSON so optional:true does not emit it to the report
        [ -s "$output_json" ] || rm -f "$output_json"
        echo "[FADE_BATCHED] Completed ${gene_id} (${direction})"
    ) &

done < "$manifest"

wait_for_all
echo "[FADE_BATCHED] Batch $batch_id finished."
