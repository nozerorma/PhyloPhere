#!/bin/bash -euo pipefail
Rscript \
'/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/subworkflows/CT/local/permulations.R' \
science.abn7829_data_s4.nex.pruned.tree \
/home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/Data/6.CAAS_traitfiles/Ranked/df8_sp.tab \
1000 \
phylogeny \
traitfile.tab \
science.abn7829_data_s4.nex.pruned.tree.resampled.output
