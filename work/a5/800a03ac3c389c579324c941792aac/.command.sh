#!/bin/bash -euo pipefail
echo "Running locally"
mkdir -p science.abn7829_data_s4.nex.pruned.resampled.output
Rscript \
'/home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/CT/local/permulations.R' \
science.abn7829_data_s4.nex.pruned.tree \
/home/miguel/IBE-UPF/PhD/NEOPLASY_PRIMATES/Malignancy_Primates/Out/2.CAAS/1.Discovery/1.Traitfiles/malignant_prevalence/malignant_prevalence.tab \
5000 \
random \
traitfile.tab \
science.abn7829_data_s4.nex.pruned.resampled.output \
500
