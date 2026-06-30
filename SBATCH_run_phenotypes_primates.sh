#!/bin/bash
#SBATCH --job-name=launch-phylophere
#SBATCH --partition=haswell
# Ensures the Slurm log directory exists
mkdir -p Slurm

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere
#SBATCH --partition=haswell
#SBATCH -t 144:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu
#SBATCH --array=1-2%2

module load Nextflow Miniconda3

source ~/.bashrc
conda deactivate || true
conda activate phylophere

# --- GLOBAL RUN CONFIG ---
export IS_TOY=false
export SOURCE_RUN_SUBDIR=""  # Set to a timestamped subdirectory name to reuse precomputed inputs, or leave empty to run from scratch
export RUN_FADE=false
export FADE_MODE="all"
export RUN_RER=false
export RER_TOOL="build_trait,build_tree,build_matrix,continuous"
export GENE_TREES="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"
export RER_PERM_BATCHES=10
export RER_PERMS_PER_BATCH=100
export RUN_VEP=true
export RUN_SCORING=true
export RUN_SCORING_STRESS=true
export RUN_CAAS_PERMULATION=true
export CAAS_FULL_PERMS=1000
export ALI_SP_NAMES="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons/results_ppga/20260615_200027/ali_sp_names.txt"
export INPUT_TAX_ID="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/5.Phylogeny/taxid_species_family_primates.tsv"

REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"
SINGLE_RUNNER="${REPO_DIR}/run_phenotype_single_primates.sh"

# --- PHENOTYPE CATALOGUE ---
case $SLURM_ARRAY_TASK_ID in
     1)  CLASS=1; TRAIT="neoplasia_prevalence"; SECONDARY="malignant_prevalence"; CTRAIT="neoplasia_necropsy"; PRUNE="neoplasia_exclude.txt"; PRUNE_SEC="malignant_exclude.txt"; DISCRETE="quintile"
         export SCORING_RER_INPUT="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/neoplasia_prevalence/runtime/filter/rerconverge/rer_results/rerconverge_summary_neoplasia_prevalence.tsv"
         export SCORING_FADE_SUMMARY_TOP="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/neoplasia_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
         export SCORING_FADE_SUMMARY_BOTTOM="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/neoplasia_prevalence/runtime/filter/selection/fade/bottom/fade_summary_bottom.tsv"
         ;;
     2)  CLASS=1; TRAIT="malignant_prevalence"; SECONDARY="neoplasia_prevalence"; CTRAIT="malignant_count";    PRUNE="malignant_exclude.txt"; PRUNE_SEC="neoplasia_exclude.txt"; DISCRETE="quintile"
         export SCORING_RER_INPUT="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/malignant_prevalence/runtime/filter/rerconverge/rer_results/rerconverge_summary_malignant_prevalence.tsv"
         export SCORING_FADE_SUMMARY_TOP="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/malignant_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
         export SCORING_FADE_SUMMARY_BOTTOM="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/diet/malignant_prevalence/runtime/filter/selection/fade/bottom/fade_summary_bottom.tsv"
         ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS $CLASS]"
echo "======================================================"

bash "$SINGLE_RUNNER" "$CLASS" "$TRAIT" "$SECONDARY" "$CTRAIT" "$PRUNE" "$PRUNE_SEC" "$DISCRETE"
EOF
}

# Submit the job
array_job_id=$(submit_array_job)
echo "Submitted array job  : ${array_job_id}  (2 tasks)"
echo "Monitor with: squeue -u \$USER"
echo "Logs in     : Slurm/"
