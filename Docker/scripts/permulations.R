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

# Set of RERconverge functions
source("simpermvec.R")

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
    # Select potential foreground and background species
    potential.fg.df <- subset(x, value %in% foreground.values)
    potential.bg.df <- subset(x, value %in% background.values)
    
    potential.fg <- potential.fg.df$name
    potential.bg <- potential.bg.df$name

    fg.species <- sample(potential.fg, foreground.size)
    bg.species <- sample(potential.bg, background.size)

  }

  else if (selection.strategy == "phylogeny") {
    # Establish tops and bottoms
    x <- x %>%
      mutate(global_label = case_when(
        value <= quantile(value, 0.25, na.rm=TRUE) ~ "low_extreme",
        value >= quantile(value, 0.75, na.rm=TRUE) ~ "high_extreme",
        TRUE ~ "normal"
      )) %>%
      ungroup()
    
    species_list_high <- x %>%
      filter(global_label == "high_extreme") %>%
      pull(name)
    species_list_low <- x %>%
      filter(global_label == "low_extreme") %>%
      pull(name)
    
    # Calculate patristic distances
    sp_distances_high <- calculate_patristic_distances(pruned.tree.o, species_list_high)
    sp_distances_low <- calculate_patristic_distances(pruned.tree.o, species_list_low)
    
    # Correct the distances
    correct_distance <- function(distance_matrix, traits, correction_type = "high") {
      species <- rownames(distance_matrix)
      
      for (i in seq_along(species)) {
        for (j in seq_along(species)) {
          if (i != j) {
            species1 <- species[i]
            species2 <- species[j]
            
            trait1 <- traits[traits$name == species1, "value"]
            trait2 <- traits[traits$name == species2, "value"]
            total <- trait1 + trait2

            # If neoplasia prevalence is 0 for both species, set the distance to its standard value squared and normalized
            if (total == 0) {
              distance_matrix[i, j] <- distance_matrix[i, j]^2
            } else {
              
              # Calculate weights
              weight1 <- trait1
              weight2 <- trait2
              
              # Apply the correction formula with weighted distances based on the correction type
              if (correction_type == "high") {
                weighted_distance <- (weight1 + weight2) * distance_matrix[i, j]^2
              } else if (correction_type == "low") {
                weighted_distance <- (2 - (weight1 + weight2)) * distance_matrix[i, j]^2
              } else {
                stop("Invalid correction type. Please choose 'high' or 'low'.")
              }
              
              # Store the weighted distance in the matrix
              distance_matrix[i, j] <- weighted_distance
            }
          }
        }
      }
      return(distance_matrix)
    }
    
    # Apply the correction function to the distance matrix
    corrected_distances_high <- correct_distance(sp_distances_high, x, correction_type = "high")
    corrected_distances_low <- correct_distance(sp_distances_low, x, correction_type = "low")
    
    # Function to rank species based on their divergence from the provided list
    rank_species <- function(distances_df, species_list) {
      ranked_species <- c()
      # Extract the most diverged species and use it as startpoint for calculating the others
      # Calculate the sum of distances for each species
      total_distances <- rowSums(distances_df)
      # Find the species with the maximum sum of distances
      max_distance_species <- names(total_distances)[which.max(total_distances)]
      # Add the species with the maximum sum of distances to the ranked list
      ranked_species <- c(ranked_species, max_distance_species)
      # Now extract rest of the species
      while (length(ranked_species) < length(species_list)) {
        temp_df <- distances_df %>%
          filter(!(rownames(.) %in% ranked_species)) %>%
          dplyr::select(all_of(ranked_species))
        # Calculate the sum of distances for each species
        total_distances <- rowSums(temp_df)
        # Find the species with the maximum sum of distances
        max_distance_species <- names(total_distances)[which.max(total_distances)]
        # Add the species with the maximum sum of distances to the ranked list
        ranked_species <- c(ranked_species, max_distance_species)
      }
      
      # Convert in ranked dataframe
      ranked_species <- data.frame(ranked_species)
      # Add numeric rank as column
      ranked_species$rank <- 1:length(ranked_species$ranked_species)
      # Rename column to species for merge
      ranked_species <- ranked_species %>%
        dplyr::rename(species = ranked_species) %>%
        # Add trait data
        left_join(x, by = c("species" = "name")) %>%
        dplyr::rename("name" = "species")
      return(ranked_species)
    }
    
    # Rank using corrected
    rankedList_high_corr <- rank_species(corrected_distances_high, species_list_high)
    rankedList_low_corr <- rank_species(corrected_distances_low, species_list_low)

    # Subset the ranked list to foreground and background species
    fg.species <- rankedList_high_corr$name[1:foreground.size]
    bg.species <- rankedList_low_corr$name[1:background.size]
  }
  
  else {

    paste("Wrong species selection options for permulations. Please select between 'random', 'inner' and 'edges'.")
  
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

