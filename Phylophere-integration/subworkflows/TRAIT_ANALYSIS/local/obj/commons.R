# Common objects for the trait analysis

## Set up variable to control command line arguments
args <- commandArgs(TRUE)


# Maybe it'd be wise to use this file to specify non-standardized attributes, such as clade or non-standard column namings
# e.g
clade_name <- "primates" # Clade of interest (used for titles, file names, etc)
taxon_of_interest <- "family" # Taxonomic level of interest (column name in the trait file)
trait <- "activity_pattern" # Trait of interest (column name in the trait file)

###################################################################################################
####### If there is a sample size trait (e.g. counts in a prevalence trait (cases/sample)), #######
####### specify it here, as it will allow for a more in-depth analysis.                     ####### 
###################################################################################################

N_trait <- "pool"

# Trait objects
## Trait data needs to have standardized colnames (species, family) and be comma-separated (csv format).
trait_df <- read.csv(args[1], sep = ",")
traits[] <- lapply(trait_df, function(col) trimws(col)) ## Remove leading/tailing whitespaces (safety)


# Phylo objects
tree <- read.newick(args[2]) # Tree object

tree_species <- tree$tip.label # Tree species

common_names <- intersect(tree$tip.label, trait_df$species) ## Tree curation
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, common_names))

# Set work directories
resultsDir <- args[3]

# Set seed for reproducibility
# Parameterize seed or set as seed local time
set.seed(args[4]) # If local time, take it from nextflow run, in any case should be parameterized
