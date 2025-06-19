#!/usr/bin/env Rscript
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
# File: permulations.R
#
# This script is part of CAASTOOLS.

# A Convergent Amino Acid Substitution identification
# and analysis toolbox
#
# Author:         Fabio Barteri (fabio.barteri@upf.edu)
#
# Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
#                 Xavier Farré (xfarrer@igtp.cat),
#                 David de Juan (david.juan@upf.edu),
#                 Miguel Ramon (miguel.ramon@upf.edu).
#
# SCRIPT NAME: permulations.r
# DESCRIPTION: Permulation script from RERconverge
# DEPENDENCIES: modules in modules/simulation folder
# CALLED BY: simulatte.py


library(tibble)
library(readr)
library(ape)
library(geiger)
library(dplyr)

# Set of RERConverge functions used by CT Resample
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm = ratematrix(treewithbranchlengths, namedvec)
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  vec = simulatedvec
  vec
}
## Simpermvec
simpermvec <- function(namedvec, treewithbranchlengths) {
  vec = simulatevec(namedvec, treewithbranchlengths)
  simsorted = sort(vec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}

# Inputs

args = commandArgs(trailingOnly=TRUE)

tree <- args[1]
config.file <- args[2]
number.of.cycles <- args[3]
selection.strategy <- args[4]
phenotypes <- args[5]
outfile <- args[6]

# tree <- "Data/5.Phylogeny/science.abn7829_data_s4.nex.pruned.tree"
# config.file <-"Out/2.CAAS/20240703-def/4.Traitfiles/df4_sp.tab"
# number.of.cycles <- "10"
# selection.strategy <- "phylogeny"
# phenotypes <- "Data/9.CAAS_bootstrap/traitfile.tab"
# outfile <- "Data/9.CAAS_bootstrap/permulated_traits.tab"

# Read the tree object.
#imported.tree <- read_file(tree)

tree.o <- read.tree(tree)
trait <- tree.o$tip.label
l <- length(trait)

# Read the config file

cfg <- read.table(config.file, sep ="\t", header = F)

foreground.df <- subset(cfg, cfg$V2 == "1")
background.df <- subset(cfg, cfg$V2 == "0")

foreground.size <- length(foreground.df$V1)
background.size <- length(background.df$V1)

foreground.species <- foreground.df$V1
background.species <- background.df$V1

# Read the phenotype file

phenotype.df <- read.table(phenotypes, sep = "\t")
foreground.values <- subset(phenotype.df, phenotype.df$V1 %in% foreground.species)$V2
background.values <- subset(phenotype.df, phenotype.df$V1 %in% background.species)$V2


# SEED AND PRUNE
starting.values <- phenotype.df$V2
all.species <- phenotype.df$V1

pruned.tree.o <- drop.tip(tree.o, setdiff(tree.o$tip.label, all.species))

names(starting.values) <- all.species

### SIMULATE
counter = 0

simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))

calculate_patristic_distances <- function(tree, species_list) {

  # Prune the tree
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, species_list))

  # Calculate patristic distances
  patristic_distances <- cophenetic(pruned_tree)

  # Convert distances to dataframe
  distances_df <- as.data.frame(as.matrix(patristic_distances))
  colnames(distances_df) <- rownames(distances_df)

  # Filter out species not in the list (if needed)
  distances_df <- distances_df[species_list, species_list]

  return(distances_df)
}


for (j in 1:as.integer(number.of.cycles)){
  counter = counter + 1
  cycle.tag = paste("b", as.character(counter), sep = "_")
  permulated_phenotype <- simpermvec(starting.values, pruned.tree.o)
  x <- enframe(permulated_phenotype)

  if (selection.strategy == "random") {
    print("Using strategy: Random")
    # Select potential foreground and background species
    potential.fg.df <- subset(x, value %in% foreground.values)
    potential.bg.df <- subset(x, value %in% background.values)

    potential.fg <- potential.fg.df$name
    potential.bg <- potential.bg.df$name

    fg.species <- sample(potential.fg, foreground.size)
    bg.species <- sample(potential.bg, background.size)

  }

  else if (selection.strategy == "phylogeny") {
    print("Using strategy: Phylogeny")
    # Establish tops and bottoms
    x <- x %>%
      dplyr::mutate(trait_label = dplyr::case_when(
        value <= quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme",
        value >= quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
        TRUE ~ "normal"
      )) %>%
      dplyr::ungroup()

    # Keep record of the species in the high and low extremes of the trait distribution
    species_list_high <- x %>%
      dplyr::filter(trait_label == "high_extreme") %>%
      dplyr::pull(name)
    species_list_low <- x %>%
      dplyr::filter(trait_label == "low_extreme") %>%
      dplyr::pull(name)
    # Combine them to work with a rectangular distance matrix
    species_list_combined <- c(species_list_high, species_list_low)

    # Calculate patristic distances
    # sp_distances_high <- calculate_patristic_distances(pruned.tree.o, species_list_high)
    # sp_distances_low <- calculate_patristic_distances(pruned.tree.o, species_list_low)
    sp_distances <- calculate_patristic_distances(pruned.tree.o, species_list_combined)

    # Correct the distances
    pair_selection <- function(distance_matrix, traits_df) {
      
      # Transform the distance matrix into a data frame
      distance_df <- as.data.frame(as.matrix(distance_matrix)) %>%
        # Separate matrix components in a long format
        rownames_to_column(var = "species1") %>%
        gather(key = "species2", value = "distance", -species1) %>%
        # Remove self-distances
        filter(distance != 0)
        # # Add trait data. Not needed now but will be in future iterations of the algorithm.
        # %>%
        # left_join(phenotype.df, by = c("species1" = "V1")) %>%
        # dplyr::rename(value1 = value) %>%
        # left_join(phenotype.df, by = c("species2" = "V1")) %>%
        # dplyr::rename(value2 = value) %>%
        ## For phenotypes with 0 values, add a small value to avoid division by 0
        # mutate(value2 = ifelse(value2 == 0, 1e-10, value2),
        ## Calculate the distance between the two species
        #   diff = abs(value1 - value2),
        ## Correct the distance using the trait value (i.e.corr_dist = diff * distance)
        ## This is naive. If phenotype is correction factor, it must be weighted to understand hwo relevant it might be.
        #   corr_dist = diff * distance)
      
      # Filter out low and high species and arrange by corrected distance
      distance_df <- distance_df %>%
        filter(!(species1 %in% species_list_low | species2 %in% species_list_high)) %>%
        arrange(distance)
      
      # Initialize results variables
      pair_reference <- data.frame()
      selected_species_top <- c()
      selected_species_bottom <- c()
      
      # REFERENCE PAIR. Find the top pair based (solely, for now) on distance
      top_pair <- distance_df %>%
        slice_head(n = 1)
      # Add top pair to the reference dataframe
      pair_reference <- rbind(pair_reference, top_pair)
      
      # Update the list of selected species
      selected_species_top <- top_pair$species1
      selected_species_bottom <- top_pair$species2
      
      # Remove rows involving selected species
      distance_df <- distance_df %>%
        filter(!(species1 %in% c(selected_species_top) | species2 %in% c(selected_species_bottom)))
      
      # Establish the cluster matrix for calculating the Dunn index
      mat <- as.matrix(distance_matrix)
      dist_mat <- as.dist(mat)
      
      # Initialize results for Dunn index
      dunn_results <- data.frame(species1 = character(), species2 = character(), Dunn_index = numeric())
      dunn_result_cummulative <- data.frame()
      
      # Start loop to find the most convergent-divergent pairs
      while (nrow(distance_df) > 0) {
        # Apply across rows to compute Dunn index for all pairs
        dunn_result <- apply(distance_df, 1, function(row) {
          query_species1 <- row["species1"]
          query_species2 <- row["species2"]
          
          # Restrict the matrix to queried species
          queried_species <- unique(c(query_species1, query_species2, selected_species_top, selected_species_bottom))
          restricted_matrix <- mat[queried_species, queried_species, drop = FALSE]
          
          # Convert the restricted matrix into a distance object
          restricted_dist <- as.dist(restricted_matrix)
          
          # Dynamically assign clusters
          cluster <- rep(2, length(queried_species))  # Default to cluster 2
          names(cluster) <- queried_species
          cluster[c(query_species1, query_species2)] <- 1
          
          # Compute Dunn index
          dunn_value <- dunn(restricted_dist, cluster, method = "complete")
          
          return(data.frame(species1 = query_species1, species2 = query_species2, Dunn_index = dunn_value))
        })
        
        # Combine results into a dataframe
        dunn_result <- do.call(rbind, dunn_result)
        
        ## DEBUG Combine into cummulative dataframe
        #dunn_result_cummulative <- rbind(dunn_result_cummulative, dunn_result)
        
        # Select the top result by Dunn index
        top_result <- dunn_result %>%
          arrange(desc(Dunn_index)) %>%
          slice_head(n = 1)
        
        # Add to Dunn results
        dunn_results <- rbind(dunn_results, top_result)
        
        # Update selected species lists
        selected_species_top <- unique(c(selected_species_top, top_result$species1))
        selected_species_bottom <- unique(c(selected_species_bottom, top_result$species2))
        
        # Remove rows involving the selected species
        distance_df <- distance_df %>%
          filter(!(species1 %in% selected_species_top | species2 %in% selected_species_bottom))
      }
      
      # Prepare output with final reference and Dunn results
      pair_reference <- pair_reference %>%
        dplyr::select(species1, species2)
      
      dunn_out <- dunn_results %>%
        dplyr::select(species1, species2)
      
      pair_out <- rbind(pair_reference, dunn_out)
      
      # Remove ugly rownames
      rownames(pair_out) <- NULL
  
      ## DEBUG 
      #return(list(dunn_results = dunn_results, pair_reference = pair_reference, dunn_result_cummulative = dunn_result_cummulative, pair_out = pair_out))
      
      return(pair_out)
    }

    # Apply the correction function to the distance matrix
    selected_pairs <- correct_distance(sp_distances, x)

    # Subset the ranked list to foreground and background species
    fg.species <- selected_pairs$species1[1:foreground.size]
    bg.species <- selected_pairs$species2[1:background.size]
  }

  else {

    paste("Wrong species selection options for permulations. Please select between 'random' and 'phylogeny'.")

  }

  # Create the output
  fg.species.tag = paste(fg.species, collapse=",")
  bg.species.tag = paste(bg.species, collapse=",")

  logline = paste("Processing permulation", as.character(counter), "of", as.character(number.of.cycles))
  outline = c(cycle.tag, fg.species.tag, bg.species.tag)
  #simulated.traits.df <- rbind(simulated.traits.df, outline)
  simulated.traits.df[nrow(simulated.traits.df) +1,] <- outline

  write(logline, stdout())
}

# Export the result
write.table(simulated.traits.df, sep="\t",col.names=FALSE, row.names = FALSE, file=outfile, quote=FALSE)

