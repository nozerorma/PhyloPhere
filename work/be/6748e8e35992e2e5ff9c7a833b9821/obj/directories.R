# ----------------------------------------
# Directory Definitions
# ----------------------------------------

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

print(paste("Working Directory:", workingDir))
print(paste("Results Directory:", resultsDir))
debug_log("Working Directory: %s", workingDir)
debug_log("Results Directory: %s", resultsDir)

# ----------------------------------------
# Common Subdirectory Paths
# ----------------------------------------

# Trait and CI directories


# Results subdirectories
data_exploration_dir <- file.path(resultsDir, "1.Data-exploration")
species_distribution_dir <- file.path(data_exploration_dir, "1.Species_distribution")
extreme_plots_dir <- file.path(data_exploration_dir, "2.Extreme_plots")
phylo_distribution_dir <- file.path(data_exploration_dir, "3.Phylogenetic_distribution")
ci_dir <- file.path(data_exploration_dir, "4.CI_overlaps")

# CAAS directories
caas_dir <- file.path(resultsDir, "2.CAAS")
traitfile_dir <- file.path(caas_dir, "1.Traitfiles")
bootstrap_traitfile_dir <- file.path(caas_dir, "2.Bootstrap_traitfiles")

debug_log("data_exploration_dir = %s", data_exploration_dir)
debug_log("species_distribution_dir = %s", species_distribution_dir)
debug_log("extreme_plots_dir = %s", extreme_plots_dir)
debug_log("phylo_distribution_dir = %s", phylo_distribution_dir)
debug_log("ci_dir = %s", ci_dir)

