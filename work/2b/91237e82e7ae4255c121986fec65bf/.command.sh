#!/bin/bash -euo pipefail
cp -R /home/miguel/IBE-UPF/PhD/PhyloPhere/Phylophere-integration/subworkflows/TRAIT_ANALYSIS/local/* .
mkdir -p reporting
Rscript -e "rmarkdown::render('1.Dataset_exploration.Rmd', output_dir='.', quiet=TRUE)" --args           'cancer_traits-FINAL-DATASET-120424.csv'           'science.abn7829_data_s4.nex.tree'           'reporting'           '1998'           'primates'           'family'           'neoplasia_prevalence'           ''           '/home/miguel/IBE-UPF/PhD/PhyloPhere/Data/5.Phylogeny/taxid_species_family.tsv'
