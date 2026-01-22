# ----------------------------------------
# Directory Definitions
# ----------------------------------------

print(paste("Working Directory:", workingDir))
print(paste("Results Directory:", resultsDir))

# ----------------------------------------
# Common Subdirectory Paths
# ----------------------------------------

# Trait and CI directories


# Results subdirectories
data_exploration_dir <- file.path(resultsDir, "1.Data-exploration")
species_distribution_dir <- file.path(data_exploration_dir, "1.Species_distribution")
extreme_plots_dir <- file.path(data_exploration_dir, "2.Extreme_plots")
phylo_distribution_dir <- file.path(data_exploration_dir, "3.Phylogenetic_distribution")
ci_overlaps_dir <- file.path(data_exploration_dir, "4.CI_overlaps")
curation_dir <- file.path(data_exploration_dir, "5.Data-curation")
post_curation_dir <- file.path(data_exploration_dir, "6.Post-curation")

# CAAS directories
caas_dir <- file.path(resultsDir, "2.CAAS")
traitfile_dir <- file.path(caas_dir, "1.Traitfiles")
bootstrap_traitfile_dir <- file.path(caas_dir, "1.5.Bootstrap_traitfiles")