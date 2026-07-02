#!/usr/bin/env bash
# run_annotate_tree_fg_batch.sh
# ───────────────────────────────────────────────────────────────────────────
# Annotate a batch of trees with foreground species labels using bash job control.
# Each gene in the manifest is launched as a background subshell; up to
# --workers run concurrently. Individual gene failures are TOLERATED (logged
# and skipped) so the overall batch task always succeeds.
#
# Manifest format (tab-separated, one gene per line):
#   gene_id <TAB> direction <TAB> fasta_filename <TAB> tree_filename <TAB> species_file_filename
#
# Files are staged by Nextflow into:
#   fastas/<fasta_filename>
#   trees/<tree_filename>
#   species_files/<species_file_filename>
# ───────────────────────────────────────────────────────────────────────────
set -uo pipefail  # NOT -e: individual gene failures must not abort the batch

batch_id=""
manifest=""
workers="1"
runner_mode=""
local_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)    batch_id="$2";    shift 2 ;;
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --runner-mode) runner_mode="$2"; shift 2 ;;
        --local-dir)   local_dir="$2";   shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$runner_mode" || -z "$local_dir" ]]; then
    echo "Missing required arguments for annotate-tree-fg batch runner" >&2
    exit 1
fi

if ! [[ "$workers" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid --workers value: $workers" >&2
    exit 1
fi

# ── Python / entrypoint ──────────────────────────────────────────────────────
if [[ "$runner_mode" == "container" ]]; then
    PY_BIN="/usr/local/bin/_entrypoint.sh python"
else
    PY_BIN="python"
fi

ANNOTATE_TREE="${local_dir}/annotate_tree_fg.py"

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "[ANNOTATE_TREE_FG_BATCHED] Batch ${batch_id}: ${gene_count} genes, ${workers} workers"

# ── Job-control helpers ──────────────────────────────────────────────────────
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
while IFS=$'\t' read -r gene_id direction fasta_name tree_name species_file_name; do
    [[ -z "${gene_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[ANNOTATE_TREE_FG_BATCHED] Launching ${gene_id} (${direction}) ($idx/$gene_count)"

    fasta_path="fastas/${fasta_name}"
    tree_path="trees/${tree_name}"
    species_file_path="species_files/${species_file_name}"

    out_tree="${gene_id}_${direction}_fg.nwk"
    out_fasta="${gene_id}_${direction}.fa"

    (
        ${PY_BIN} "${ANNOTATE_TREE}" \
            --species-file "${species_file_path}" \
            --tree         "${tree_path}" \
            --fasta        "${fasta_path}" \
            --fasta_out    "${out_fasta}" \
            --output       "${out_tree}" \
            || echo "[ANNOTATE_TREE_FG_BATCHED] annotate_tree_fg failed for ${gene_id} (${direction}), skipping"

        # Remove empty outputs so optional:true does not emit them
        [ -s "${out_tree}" ] || rm -f "${out_tree}"
        [ -s "${out_fasta}" ] || rm -f "${out_fasta}"

        echo "[ANNOTATE_TREE_FG_BATCHED] Completed ${gene_id} (${direction})"
    ) &

done < "$manifest"

wait_for_all
echo "[ANNOTATE_TREE_FG_BATCHED] Batch $batch_id finished."
