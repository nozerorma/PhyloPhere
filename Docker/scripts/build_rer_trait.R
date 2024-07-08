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
# File: build_rer_trait.R
#

# Set up variable to control command line arguments
args <- commandArgs(TRUE)

# Load libraries
library(dplyr)
library(RERconverge)

# Load traitfile
## Pass as arg from nextflow config, should be
cancer_traits <- read.csv(args[1])

## Remove whitespaces
cancer_traits[] <- lapply(cancer_traits, function(col) trimws(col))

## Binarize names
cancer_traits$species <- gsub(" ", "_", cancer_traits$species)

## Reorder and filter species
cancer_traits <- cancer_traits %>%
  select(species,neoplasia_prevalence) %>%
  filter(!is.na(neoplasia_prevalence))

# Select only those columns of interest and build named vector
neoplasia_vector <- setNames(as.numeric(cancer_traits$neoplasia_prevalence), cancer_traits$species)
head(neoplasia_vector)

# Write the generated vector
## Parameterize to traits_continuous
neoplasiaPath <- args[2]
print(neoplasiaPath)
save(neoplasia_vector, file = neoplasiaPath)

### DONE ###
