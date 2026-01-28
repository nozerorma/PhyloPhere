#!/bin/bash -euo pipefail
cp -R /home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/TRAIT_ANALYSIS/local/* .
mkdir -p reporting/HTML_reports
Rscript -e "rmarkdown::render('2.Phenotype_exploration.Rmd', output_dir='reporting/HTML_reports', quiet=TRUE)" --args           'cancer_traits_processed-LQ.csv'           'science.abn7829_data_s4.nex.tree'           'reporting'           '1998'           'primates'           'family'           'neoplasia_prevalence'           'adult_necropsy_count'           'neoplasia_necropsy'           '/home/miguel/IBE-UPF/PhD/PhyloPhere/Data/5.Phylogeny/taxid_species_family.tsv'           'malignant_prevalence'           'LQ'
