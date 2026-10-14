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
# Manifest (tab-separated, one run per line, written by the process):
#   gene_id <TAB> direction <TAB> fasta_filename <TAB> annotated_tree_filename
# A 3-column line (gene_id, fasta, tree) takes its direction from --direction.
# The files are staged by Nextflow as fastas/<fasta_filename> and trees/<annotated_tree_filename>.
# The runs of one gene share the baseline model fit (--share-baseline true) through a per-gene cache file.
#
# Args (named flags):
#   --batch-id, --manifest, --direction (top | bottom | both), --runner-mode (container | local)
#       required; container mode runs hyphy through /usr/local/bin/_entrypoint.sh
#   --workers            concurrent genes (positive integer, default 1)
#   --cpu-per-worker     CPUs of each HyPhy call (default 1)
#   --model, --model-file-arg, --method, --grid, --concentration
#       FADE options (--model-file-arg is split into separate tokens)
#   --mcmc-chains, --mcmc-chain-length, --mcmc-burn-in, --mcmc-samples
#       MCMC options, used only when --method is not Variational-Bayes
#   --bf-threshold      Bayes factor threshold of a significant site (default 100)
# Each HyPhy JSON is reduced to a gene row and the significant sites of its gene by fade_json_to_tables.py and then deleted; a gene
# without a usable result writes no row.
# Output: <batch_id>.fade_genes.tsv and <batch_id>.fade_sites.tsv in the working directory (header only when no gene has a result).
# =============================================================================
set -uo pipefail  # no -e: a failed gene must not abort the batch

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fade_bf="${script_dir}/FADE_phylophere.bf"
extract_py="${script_dir}/fade_json_to_tables.py"

batch_id=""
manifest=""
direction="both"
share_baseline="true"
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
bf_threshold="100"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)           batch_id="$2";          shift 2 ;;
        --manifest)           manifest="$2";           shift 2 ;;
        --direction)          direction="$2";          shift 2 ;;
        --share-baseline)     share_baseline="$2";     shift 2 ;;
        --fade-bf)            fade_bf="$2";            shift 2 ;;
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
        --bf-threshold)       bf_threshold="$2";       shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$runner_mode" ]]; then
    echo "Missing required arguments for FADE batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value: $workers" >&2
    exit 1
fi

# ── Base HyPhy command ───────────────────────────────────────────────────────
declare -a base_cmd
declare -a py_cmd
if [[ "$runner_mode" == "container" ]]; then
    base_cmd=("/usr/local/bin/_entrypoint.sh" "hyphy" "$fade_bf")
    py_cmd=("/usr/local/bin/_entrypoint.sh" "python")
else
    base_cmd=("hyphy" "$fade_bf")
    py_cmd=("python")
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

row_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "Running batched FADE task $batch_id (direction: $direction, share_baseline: $share_baseline)"
echo "Runs in batch: $row_count"
echo "Concurrent workers: $workers (${cpu_per_worker} CPU(s) per HyPhy call)"
mkdir -p runs   # one pair of tables per run, concatenated at the end of the batch

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

launch_gene_group() {
    local gid="$1"
    shift
    local -a specs=("$@")
    [[ "${#specs[@]}" -eq 0 ]] && return 0

    wait_for_slot
    (
        export OMP_NUM_THREADS="${cpu_per_worker}"
        export MKL_NUM_THREADS="${cpu_per_worker}"
        export OPENBLAS_NUM_THREADS="${cpu_per_worker}"
        export BLAS_NUM_THREADS="${cpu_per_worker}"

        local cache_arg="NONE"
        if [[ "${#specs[@]}" -gt 1 && "$share_baseline" == "true" ]]; then
            # Absolute path: HyPhy resolves a relative cache path against the directory of the batch file,
            # not the working directory, so a relative name would be written beside FADE_phylophere.bf.
            cache_arg="${PWD}/${gid}.baseline.cache"
            rm -f "$cache_arg"
        fi
        trap 'rm -f "${gid}.baseline.cache" "${gid}".*.hyphy.log' EXIT

        for spec in "${specs[@]}"; do
            IFS='|' read -r run_dir fasta_name tree_name <<< "$spec"
            local output_json="${gid}.${run_dir}.FADE.json"
            local fasta_path="fastas/${fasta_name}"
            local tree_path="trees/${tree_name}"
            local tmp_log="${gid}.${run_dir}.hyphy.log"

            echo "[FADE_BATCHED] Launching ${gid} (${run_dir})"

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
                "--cache"                   "$cache_arg"
                "--output"                  "$output_json"
            )
            if [[ -n "$model_file_arg" ]]; then
                read -r -a _mfa <<< "$model_file_arg"
                cmd+=("${_mfa[@]}")
            fi
            cmd+=("${mcmc_args[@]}")

            if ! "${cmd[@]}" > "$tmp_log" 2>&1; then
                echo "[FADE_BATCHED] FADE failed for ${gid} (${run_dir}), skipping"
                tail -n 25 "$tmp_log" >&2 || true
            fi
            # Marker kept in the task log (the HyPhy log is discarded): the shared baseline fit was reused.
            if grep -q "Loaded baseline model fit from cache" "$tmp_log" 2>/dev/null; then
                echo "[FADE_BATCHED] ${gid} (${run_dir}): baseline model fit loaded from the shared cache"
            fi
            rm -f "$tmp_log" "${fasta_path}.FADE.cache"

            # The tables of the run replace the JSON, which is not kept: a failed or empty run writes no row.
            if [[ -s "$output_json" ]]; then
                "${py_cmd[@]}" "$extract_py" --json "$output_json" --gene "$gid" --direction "$run_dir" \
                    --bf-threshold "$bf_threshold" \
                    --genes-out "runs/${gid}.${run_dir}.genes.tsv" --sites-out "runs/${gid}.${run_dir}.sites.tsv" \
                    || echo "[FADE_BATCHED] ${gid} (${run_dir}): no usable FADE result, skipping" >&2
            fi
            rm -f "$output_json"
            echo "[FADE_BATCHED] Completed ${gid} (${run_dir})"
        done

        rm -f "${gid}.baseline.cache"
    ) &
}

# ── Group consecutive manifest rows by gene_id and launch ────────────────────
current_gene=""
declare -a current_specs=()

while IFS=$'\t' read -r col1 col2 col3 col4; do
    [[ -z "${col1:-}" ]] && continue
    if [[ -n "${col4:-}" ]]; then
        gene_id="$col1"
        run_dir="$col2"
        fasta_name="$col3"
        tree_name="$col4"
    else
        gene_id="$col1"
        run_dir="$direction"
        fasta_name="$col2"
        tree_name="$col3"
    fi

    if [[ -n "$current_gene" && "$gene_id" != "$current_gene" ]]; then
        launch_gene_group "$current_gene" "${current_specs[@]}"
        current_specs=()
    fi
    current_gene="$gene_id"
    current_specs+=("${run_dir}|${fasta_name}|${tree_name}")
done < "$manifest"

if [[ -n "$current_gene" && "${#current_specs[@]}" -gt 0 ]]; then
    launch_gene_group "$current_gene" "${current_specs[@]}"
fi

wait_for_all
rm -f ./*.baseline.cache fastas/*.FADE.cache

# One table of each kind per batch: the header once, then the rows of every run.
for kind in genes sites; do
    out="${batch_id}.fade_${kind}.tsv"
    "${py_cmd[@]}" "$extract_py" --print-header "$kind" > "$out"
    for f in runs/*."${kind}".tsv; do
        [[ -e "$f" ]] && tail -n +2 "$f" >> "$out"
    done
done
rm -rf runs
echo "[FADE_BATCHED] Batch $batch_id finished."
