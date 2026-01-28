#!/bin/bash -euo pipefail
cp -R /home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/TRAIT_ANALYSIS/local/* .
mkdir -p data_exploration/HTML_reports
Rscript -e "rmarkdown::render('4.Independent_contrasts.Rmd', output_dir='data_exploration/HTML_reports', quiet=TRUE)" --args           'cancer_traits_processed-LQ.csv'           'science.abn7829_data_s4.nex.tree'           'data_exploration'           '1998'           'primates'           'family'           'neoplasia_prevalence'           'adult_necropsy_count'            'neoplasia_necropsy'           '/home/miguel/IBE-UPF/PhD/PhyloPhere/Data/5.Phylogeny/taxid_species_family.tsv'           'malignant_prevalence'           'LQ'           ''
