#!/usr/bin/env bash
# run_prep_alignments_batch.sh
# ───────────────────────────────────────────────────────────────────────────
# Convert and tree-filter a batch of alignment files in a single Nextflow task.
# For each gene in the manifest:
#   1. Detect format (FASTA or PHYLIP-like)
#   2. Convert to FASTA if needed (phylip_to_fasta.py)
#   3. Filter sequences to tree taxa (filter_fasta_to_tree.py)
#   4. Write <gene_id>.filtered.fa
#
# Individual gene failures are TOLERATED (logged and skipped) so the overall
# batch task always succeeds.
#
# Manifest format (tab-separated, one gene per line):
#   gene_id <TAB> alignment_filename
#
# Alignment files are staged by Nextflow under:  alignments/<alignment_filename>
# ───────────────────────────────────────────────────────────────────────────
set -uo pipefail  # NOT -e: individual gene failures must not abort the batch

batch_id=""
manifest=""
workers="4"
runner_mode=""
tree=""
local_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch-id)    batch_id="$2";   shift 2 ;;
        --manifest)    manifest="$2";   shift 2 ;;
        --workers)     workers="$2";    shift 2 ;;
        --runner-mode) runner_mode="$2"; shift 2 ;;
        --tree)        tree="$2";       shift 2 ;;
        --local-dir)   local_dir="$2";  shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$batch_id" || -z "$manifest" || -z "$runner_mode" || -z "$tree" || -z "$local_dir" ]]; then
    echo "Missing required arguments for prep-alignments batch runner" >&2
    echo "  batch_id='${batch_id}' manifest='${manifest}' runner_mode='${runner_mode}' tree='${tree}' local_dir='${local_dir}'" >&2
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

PHYLIP_TO_FASTA="${local_dir}/phylip_to_fasta.py"
FILTER_FASTA="${local_dir}/filter_fasta_to_tree.py"

gene_count="$(grep -cve '^[[:space:]]*$' "$manifest" || true)"
echo "[PREP_ALIGNMENTS_BATCHED] Batch ${batch_id}: ${gene_count} genes, ${workers} workers"

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
while IFS=$'\t' read -r gene_id alignment_name; do
    [[ -z "${gene_id:-}" ]] && continue
    idx=$((idx + 1))
    wait_for_slot
    echo "[PREP_ALIGNMENTS_BATCHED] Launching ${gene_id} (${idx}/${gene_count})"

    alignment_path="alignments/${alignment_name}"
    tmp_fasta="${gene_id}_tmp.fa"
    out_fasta="${gene_id}.filtered.fa"

    (
        # Step 1: convert to FASTA (or copy if already FASTA)
        _lower="$(echo "${alignment_name}" | tr '[:upper:]' '[:lower:]')"
        if [[ "$_lower" == *.fa || "$_lower" == *.fasta ]]; then
            cp "${alignment_path}" "${tmp_fasta}"
        else
            ${PY_BIN} "${PHYLIP_TO_FASTA}" "${alignment_path}" "${tmp_fasta}" \
                || { echo "[PREP_ALIGNMENTS_BATCHED] PHYLIP→FASTA failed for ${gene_id}, skipping" >&2; exit 0; }
        fi

        # Step 2: filter to tree taxa
        ${PY_BIN} "${FILTER_FASTA}" \
            --tree   "${tree}" \
            --fasta  "${tmp_fasta}" \
            --output "${out_fasta}" \
            || { echo "[PREP_ALIGNMENTS_BATCHED] filter_fasta_to_tree failed for ${gene_id}, skipping" >&2; rm -f "${out_fasta}"; exit 0; }

        # Remove empty outputs so optional:true does not emit them
        [ -s "${out_fasta}" ] || rm -f "${out_fasta}"

        rm -f "${tmp_fasta}"
        echo "[PREP_ALIGNMENTS_BATCHED] Completed ${gene_id}"
    ) &

done < "$manifest"

wait_for_all
echo "[PREP_ALIGNMENTS_BATCHED] Batch ${batch_id} finished."
