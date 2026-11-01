#!/bin/bash
#SBATCH --job-name=def2_pp2a
#SBATCH --partition=high-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=def2_pp2a-%j.out
#SBATCH --error=def2_pp2a-%j.err
set -Eeuo pipefail

# Submit from the synchronized PhyloPhere repository on correfoc.
REPO_DIR="${PHYLOPHERE_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
source /etc/profile
set +u
module load Miniconda3 2>/dev/null || module load Miniconda3/25.7.0-2 2>/dev/null || true
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate phylophere
fi
set -u
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
OLD_ROOT="${OLD_ROOT:-/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/malignant_def2/malignant_prevalence_complete}"
NEW_BASE="${NEW_BASE:-/homes/users/mramon/scratch/3.Work_dirs_malignant_def}"
OUT_ROOT="${OUT_ROOT:-/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/def2_pp2a_updated_U_${SLURM_JOB_ID}}"
cd "$REPO_DIR"
python3 bin/prepare_def2_pp2a.py --repo "$REPO_DIR" --old-root "$OLD_ROOT" \
    --new-work "$NEW_BASE/malignant_prevalence_complete" \
    --new-params "$NEW_BASE/malignant_prevalence_complete_nxf_run/params.json" \
    --output "$OUT_ROOT" --workers "${SLURM_CPUS_PER_TASK:-4}"

# Explicit legacy consumer, compatible with def2's stored score/block schemas.
# The new universe is already prepared. This does not bypass or forge v4 provenance.
mkdir "$OUT_ROOT/consumer_source"
git archive 5dc078ce subworkflows/POSENRICH/local/src | tar -x -C "$OUT_ROOT/consumer_source"
POS_SRC="$OUT_ROOT/consumer_source/subworkflows/POSENRICH/local/src"
shopt -s nullglob
blocks=("$OLD_ROOT"/caas_permulation/blocks/*)
if (( ${#blocks[@]} == 0 )); then echo 'No def2 null blocks found' >&2; exit 1; fi
python3 "$POS_SRC/posenrich_within_gene.py" \
    --universe "$OUT_ROOT/universe/universe.tsv.gz" --observed "$OUT_ROOT/position_scores.tsv" \
    --registry-dir "$OUT_ROOT/registry" --blocks "${blocks[@]}" \
    --collections custom --min-genes 2 --holdout 0 --output-dir "$OUT_ROOT/within_gene"
python3 "$POS_SRC/posenrich_proximity.py" \
    --registry_dir "$OUT_ROOT/registry" --output_dir "$OUT_ROOT/proximity" \
    --collections custom --radios 0,5,10
python3 - "$OUT_ROOT" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1]);audit=json.loads((root/'recovery_audit.json').read_text())
rows=list(csv.DictReader((root/'within_gene/within_gene_results.tsv').open(),delimiter='\t'))
text=['# PP2A: def2 scores and updated positional universe','',
      f"Universe: {audit['universe_positions']} positions in the def2 cleaned gene background.",
      f"Excluded old observed rows outside updated universe: {audit['observed_rows_outside_updated_universe']}.",'',
      'ASR, scores and null filters are frozen to def2. Position p-values/FDR were not reused.',
      'The p/q below are the freshly computed within-gene tests; q is adjusted over the PP2A terms tested.',
      'These phosphosite sets describe association with PP2A perturbations, not proof of direct substrates.','',
      '| Term | Eligible genes | p | q | Frequency difference |','|---|---:|---:|---:|---:|']
for r in rows:
    text.append('| '+' | '.join(r.get(k,'NA') for k in ['term_id','n_genes_eligible','p','q','frequency_difference'])+' |')
text+=['','Gene contributions: within_gene/within_gene_by_gene.tsv.',
       'Exact overlaps / near positions: proximity/proximity_protein_hits.tsv.',
       'Distances are descriptive, without an enrichment p-value.']
(root/'PP2A_summary.md').write_text('\n'.join(text)+'\n')
PY
echo "PP2A output: $OUT_ROOT"
echo 'Primary result: within_gene/within_gene_results.tsv (p/q across PP2A terms)'
echo 'Gene contribution: within_gene/within_gene_by_gene.tsv'
echo 'Proximity: proximity/proximity_protein_hits.tsv'
