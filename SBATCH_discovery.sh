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
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: SBATCH_discovery.sh
#

#SBATCH --job-name=nfct-discovery
#SBATCH -p haswell
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=8G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --time=8:00:00

# send mail if needed
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu

#Define modules
module load Nextflow

# Define the directory where trait files are located
TRAIT_DIR="/gpfs42/robbyfs/scratch/lab_anavarro/mramon/nf_caastools/Data/Traitfiles/"

# Loop through trait files in the directory with a .tab file extension
for TRAIT_FILE in "$TRAIT_DIR"*.tab
do
        # Run caastools in Nextflow using the current trait file
        srun -n1 --exclusive nextflow run main.nf -with-singularity -with-tower -profile singularity --ct_tool discovery --traitfile "$TRAIT_FILE" &
done

wait
