#!/usr/bin/env bash
set -euo pipefail

# CT disambiguation-only Nextflow launcher (stub mode by default)

ROOT_DIR="/home/miguel/IBE-UPF/PhD/PhyloPhere"
MAIN_NF="${ROOT_DIR}/main.nf"

PROFILE="${PROFILE:-local}"
OUTDIR="${OUTDIR:-${ROOT_DIR}/Out_stub_ct_disambiguation}"
WORKDIR="${WORKDIR:-${ROOT_DIR}/work_stub_ct_disambiguation}"

# Keep this script focused exclusively on CT disambiguation
nextflow run "${MAIN_NF}" \
  -profile "${PROFILE}" \
  -w "${WORKDIR}" \
  --outdir "${OUTDIR}" \
  --ct_disambiguation true \
  --caas_postproc false \
  --caas_signification false \
  --contrast_selection false \
  "$@"
