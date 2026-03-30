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

workingDir <- getwd()
print(workingDir)

# Load libraries
library(dplyr)
library(RERconverge)

# Load traitfile
## File must have a multi-column structure with at least a "species" and "trait" column
## Auto-detect separator (tab or comma) by inspecting the first line
.first_line <- readLines(args[1], n = 1)
.sep <- if (grepl("\t", .first_line)) "\t" else ","
ori_traits <- read.table(args[1], header = TRUE, sep = .sep, check.names = FALSE, quote = "\"")

## Remove whitespaces
ori_traits[] <- lapply(ori_traits, function(col) trimws(col))

# Specify column name for species
sp_colname <- args[2]

## Binarize names if not performed already
ori_traits$species <- gsub(" ", "_", ori_traits[[sp_colname]])

# Reorder and filter species
selected_column <- args[3]

ori_traits <- ori_traits %>%
  dplyr::select(species, !!as.name(selected_column)) %>%
  dplyr::filter(!is.na(!!as.name(selected_column)))

# Select only those columns of interest and build named vector
trait_vector <- setNames(as.numeric(ori_traits[[selected_column]]), ori_traits[[sp_colname]])

# Write the generated vector
## Parameterize to traits_continuous
traitPath <- args[4]
print(traitPath)
save(trait_vector, file = traitPath)

# ── Detect trait type (binary vs continuous) from value distribution ──────────
# Binary: exactly 2 unique non-NA values (expected: 0/1 encoding).
# If values are not 0/1, recode min→0 / max→1 with a warning.
unique_vals <- sort(unique(na.omit(trait_vector)))
n_unique    <- length(unique_vals)

if (n_unique == 2) {
  if (!all(unique_vals == c(0, 1))) {
    warning(sprintf(
      "[RER_TRAIT] Binary values are not 0/1 (found: %s, %s). Recoding min→0, max→1.",
      unique_vals[1], unique_vals[2]
    ))
    trait_vector <- ifelse(trait_vector == unique_vals[1], 0, 1)
    # Overwrite saved trait with recoded version
    save(trait_vector, file = traitPath)
  }
  trait_type <- "binary"
  message(sprintf(
    "[RER_TRAIT] Binary trait detected (%d foreground / %d background species).",
    sum(trait_vector == 1, na.rm = TRUE),
    sum(trait_vector == 0, na.rm = TRUE)
  ))
} else {
  trait_type <- "continuous"
  message(sprintf(
    "[RER_TRAIT] Continuous trait detected (%d unique values, range [%.4g, %.4g]).",
    n_unique, min(trait_vector, na.rm = TRUE), max(trait_vector, na.rm = TRUE)
  ))
}

# Write trait type file so Nextflow can route to the correct analysis process
typeOutPath <- args[5]
writeLines(trait_type, con = typeOutPath)
message(sprintf("[RER_TRAIT] Trait type '%s' written to: %s", trait_type, typeOutPath))

### DONE ###
