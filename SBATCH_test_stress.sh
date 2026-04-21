#!/bin/bash
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
#
# PHYLOPHERE: Stress-Test SLURM Array Controller
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Submits one SLURM array task per trait, each calling test_stress_single.sh.
# All toggle flags (integrated/standalone modes, sub-features) are set ONCE
# here and exported to the per-task runner. Trait-specific parameters are
# injected via the case statement at the bottom.
#
# ─── USAGE ────────────────────────────────────────────────────────────────────
#   bash SBATCH_test_stress.sh
#
# To run only a subset of traits, change --array=1-10 to e.g. --array=1,2.
# To run a single trait outside SLURM:
#   export TRAIT="neoplasia_prevalence" TRAIT_CLASS=1 SECONDARY_TRAIT="malignant_prevalence" ...
#   bash test_stress_single.sh
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: SBATCH_test_stress.sh
#

# Ensure the Slurm log directory exists
mkdir -p Slurm

# ============================================================
# CONFIGURATION  (edit before submitting)
# ============================================================
REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"

# ─── Phenotype catalogue (array 1-10) ─────────────────────────────────────────
#   1.  neoplasia_prevalence  [CLASS 1]
#   2.  malignant_prevalence  [CLASS 1]
#   3.  frug_idx              [CLASS 2]
#   4.  fol_idx               [CLASS 2]
#   5.  ins_idx               [CLASS 2]
#   6.  omn_idx               [CLASS 2]
#   7.  omn_spec_idx          [CLASS 2]
#   8.  herb_idx              [CLASS 2]
#   9.  folfrug_idx           [CLASS 2]
#  10.  Ethanol               [CLASS 2]

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere-stress
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

module purge; module load modulepath/haswell
module load Nextflow
module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate phylophere

# ═══════════════════════════════════════════════════════════════════════════════
#   MASTER TOGGLES  —  same for every trait in this submission
#   Export them so test_stress_single.sh picks them up via ${VAR:-default}.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Infrastructure ────────────────────────────────────────────────────────────
export CLEAN_WORK=false

# ── Integrated run toggles ────────────────────────────────────────────────────
export RUN_INT_FILTER=false
export RUN_INT_EXPLORATORY=false
export INT_RESUME=false

export INT_PRUNE_DATA=false
export INT_REPORTING=false
export INT_CONTRAST_SELECTION=false
export INT_CT_POSTPROC=true
export INT_CT_SIGNIFICATION=false
export INT_CT_DISAMBIGUATION=false
export INT_ORA=true
export INT_STRING=true
export INT_CT_ACCUMULATION=false
export INT_VEP=false
export INT_FADE=false
export FADE_MODE="gene_set"
export INT_RER=true
export INT_SCORING=true
export INT_SCORING_STRESS=true
export INT_SCORING_STRESS_TOP_N=25

export INT_USE_SECONDARY_TRAIT=true
export INT_USE_BRANCH_TRAIT=true
export INT_USE_N_TRAIT=true
export INT_USE_C_TRAIT=true

export INT_DISAMBIG_ASR_MODE="precomputed"
export INT_DISAMBIG_ASR_CACHE_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data/asr"

# ── Standalone run toggles ────────────────────────────────────────────────────
export RUN_SA_CT=false
export RUN_SA_PRUNE=false
export SA_PRUNE_ENABLE=false
export RUN_SA_SIGNIFICATION=false
export RUN_SA_DISAMBIGUATION=false
export SA_DISAMBIG_COMPUTE=false
export SA_DISAMBIG_PRECOMPUTED=true
export RUN_SA_POSTPROC_FILTER=false
export RUN_SA_POSTPROC_EXPLORATORY=false
export RUN_SA_ORA=false
export RUN_SA_STRING=false
export RUN_SA_ACCUMULATION=false
export RUN_SA_REPORTING=false
export RUN_SA_CONTRAST_SELECTION=false
export RUN_SA_VEP=false
export RUN_SA_SCORING=false
export SA_SCORING_STRESS=false
export SA_SCORING_STRESS_TOP_N=25
export RUN_SA_FADE=false
export SA_FADE_MODE="gene_set"
export RUN_SA_RER=true

export SA_RER_TOOL="build_trait,build_tree,build_matrix,continuous"
export SA_RER_CONTINUOUS_ONLY=false
export SA_RER_PERM_BATCHES=100
export SA_RER_PERMS_PER_BATCH=100
export SA_RER_PERM_MODE="cc"
export SA_RER_GMT_FILE="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere/subworkflows/RERCONVERGE/dat/c2.cp.pid.v2026.1.Hs.symbols.gmt"

# ── Toy mode ──────────────────────────────────────────────────────────────────
export TOY_MODE=false
export TOY_N=1000

# ── Cluster / environment paths (override defaults in test_stress_single.sh) ──
export DATADIR="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/1.Primates_data"
export CAAS_OUTBASE="/data/samanthafs/scratch/lab_anavarro/mramon/2.Primates/2.Primates_results/CAAS_RESULTS/final"
export WORK_BASE="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs/final"
export ALI_ARCHIVE_ROOT="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons/alignments/Primates_BMGE"
export ALI_DIR="${ALI_ARCHIVE_ROOT}/PROT"
export INPUT_GENE_TREES="${DATADIR}/3.Gene_trees/Gene_trees/ALL_FEB23_geneTrees.txt"
export INPUT_ASR_CACHED="${DATADIR}/asr"
export INPUT_VEP_CDS_DIR="${ALI_ARCHIVE_ROOT}/TRIM"
export INPUT_VEP_TRACK_DIR="${ALI_ARCHIVE_ROOT}/TRACK"
export INPUT_VEP_PRIMATEAI_DB="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere/subworkflows/VEP/dat/PrimateAI-3D.hg38.txt.gz"
export INPUT_VEP_TRANSVAR_REFERENCE="/homes/users/mramon/scratch/0.Phylophere/subworkflows/VEP/dat/transvar/hg38.fa"
# Optional: precomputed VEP stage outputs to skip aa2nuc/aa2prot.
# Leave empty to keep default behavior.
export INPUT_VEP_AA2NUC_INPUT="${INPUT_VEP_AA2NUC_INPUT:-}"
export INPUT_VEP_AA2PROT_INPUT="${INPUT_VEP_AA2PROT_INPUT:-}"

export SOURCE_RUN_SUBDIR="runtime"
export CYCLES="1000000"

REPO_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/0.Phylophere"
SINGLE_RUNNER="${REPO_DIR}/test_stress_single.sh"

# ═══════════════════════════════════════════════════════════════════════════════
#   TRAIT CATALOGUE
#   Args exported to test_stress_single.sh via environment:
#     TRAIT_CLASS  TRAIT  SECONDARY_TRAIT  N_TRAIT  C_TRAIT  BRANCH_TRAIT
#     PRUNE_FILE   PRUNE_FILE_SECONDARY
# ═══════════════════════════════════════════════════════════════════════════════
# CLASS 1 traits use the cancer CSV + prune lists.
# CLASS 2 traits use the diet CSV; secondary/n/c/prune fields are left empty.
case $SLURM_ARRAY_TASK_ID in
     1)  export TRAIT_CLASS=1
         export TRAIT="neoplasia_prevalence"
         export SECONDARY_TRAIT="malignant_prevalence"
         export N_TRAIT="adult_necropsy_count"
         export C_TRAIT="neoplasia_necropsy"
         export BRANCH_TRAIT="LQ"
         export PRUNE_FILE="neoplasia_exclude.txt"
         export PRUNE_FILE_SECONDARY="malignant_exclude.txt"
         ;;
     2)  export TRAIT_CLASS=1
         export TRAIT="malignant_prevalence"
         export SECONDARY_TRAIT="neoplasia_prevalence"
         export N_TRAIT="adult_necropsy_count"
         export C_TRAIT="malignant_count"
         export BRANCH_TRAIT="LQ"
         export PRUNE_FILE="malignant_exclude.txt"
         export PRUNE_FILE_SECONDARY="neoplasia_exclude.txt"
         ;;
#      3)  export TRAIT_CLASS=2; export TRAIT="frug_idx";     export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      4)  export TRAIT_CLASS=2; export TRAIT="fol_idx";      export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      5)  export TRAIT_CLASS=2; export TRAIT="ins_idx";      export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      6)  export TRAIT_CLASS=2; export TRAIT="omn_idx";      export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      7)  export TRAIT_CLASS=2; export TRAIT="omn_spec_idx"; export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      8)  export TRAIT_CLASS=2; export TRAIT="herb_idx";     export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#      9)  export TRAIT_CLASS=2; export TRAIT="folfrug_idx";  export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
#     10)  export TRAIT_CLASS=2; export TRAIT="Ethanol";      export SECONDARY_TRAIT=""; export N_TRAIT=""; export C_TRAIT=""; export BRANCH_TRAIT="LQ"; export PRUNE_FILE=""; export PRUNE_FILE_SECONDARY="" ;;
esac

echo "======================================================"
echo " SLURM TASK ID : $SLURM_ARRAY_TASK_ID"
echo " TRAIT         : $TRAIT  [CLASS ${TRAIT_CLASS}]"
echo " NODE          : $(hostname)"
echo "======================================================"

bash "$SINGLE_RUNNER"
EOF
}

# Main
array_job_id=$(submit_array_job)
echo "Submitted array job  : ${array_job_id}  (10 tasks)"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Logs in     : Slurm/"
