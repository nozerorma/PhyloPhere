#!/bin/bash -euo pipefail
cp -R /home/miguel/IBE-UPF/PhD/PhyloPhere/Phylophere-integration/subworkflows/TRAIT_ANALYSIS/local/* .
mkdir -p reporting
Rscript -e "rmarkdown::render('1.Dataset_exploration.Rmd', output_dir='.', quiet=TRUE)" --args           'traitfile.tab'           'tree.nwk'           'reporting'           ''           'clade'           'family'           'trait'           ''
