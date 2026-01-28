#!/bin/bash -euo pipefail
cp -R /home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/TRAIT_ANALYSIS/local/* .
mkdir -p CI_composition_report/HTML_reports
Rscript -e "rmarkdown::render('3.CI-composition.Rmd', output_dir='CI_composition_report/HTML_reports', quiet=TRUE)" --args           'cancer_traits_processed-LQ.csv'           'science.abn7829_data_s4.nex.tree'           'CI_composition_report'           '1998'           'primates'           'family'           'neoplasia_prevalence'           'adult_necropsy_count'            'neoplasia_necropsy'           '/home/miguel/IBE-UPF/PhD/PhyloPhere/Data/5.Phylogeny/taxid_species_family.tsv'           'malignant_prevalence'           'LQ'
