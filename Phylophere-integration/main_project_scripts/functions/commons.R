### Common objects for the primate neoplasia project ###
# MIGUEL RAMON ALONSO

# Data import and curation

## Tree import
# Set phylodir
phylo_dir <- file.path(dataDir, "5.Phylogeny/")
# Open tree file
primate_tree <- read.newick(file.path(phylo_dir,"science.abn7829_data_s4.nex.tree"))
# Extract tree species
tree_species <- primate_tree$tip.label
names_w_homo <- c(tree_species, "Homo_sapiens")

# Set cancer trait dir
## Folder containing trait data
trait_dir <- file.path(dataDir, "1.Cancer_data/Neoplasia_species360")

# Traits
## Data should have standardized colnames (species, family). Specific traitnames will prob
## need to be hardcoded. Little by little.
traits <- read.csv(file.path(trait_dir, "traits_with_family.csv"), sep = ",")
# Remove leading/tailing whitespaces
traits[] <- lapply(traits, function(col) trimws(col))

# Curate using tax_id
# Family data
family_data <- read_tsv(file.path(dataDir, "5.Phylogeny/taxid_species_family.tsv"))

# Your primate tree
primate_tree <- as.phylo(primate_tree)

# Map species in the tree to the tax_ids
tree_ids <- data.frame(
  species = primate_tree$tip.label,
  tax_id = sapply(primate_tree$tip.label, function(x) {
    row <- family_data %>% filter(gsub(" ", "_", species) == x)
    if (nrow(row) == 1) {
      return(row$tax_id)
    } else {
      return(NA)
    }
  })
)

# Data curation
## Curate dataset using tree species
common_tax_ids <- intersect(tree_ids$tax_id, traits$tax_id)

# Prune tree to species for which we have phenotype
pruned_tree <- drop.tip(primate_tree, primate_tree$tip.label[!tree_ids$tax_id %in% common_tax_ids])