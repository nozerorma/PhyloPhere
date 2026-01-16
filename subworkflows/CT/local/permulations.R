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
outdir <- args[6]
chunk.size <- ifelse(length(args) >= 7, as.integer(args[7]), 500)

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  write(paste("[INFO]", Sys.time(), "Created output directory:", outdir), stdout())
}

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
chunk.counter = 1
file.counter = 1

simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))

# Start timing
start.time <- Sys.time()
write(paste("[START]", start.time, "Beginning permulation generation..."), stdout())
write(paste("[INFO] Total cycles:", number.of.cycles, "| Chunk size:", chunk.size), stdout())

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
        # Round distance and diff
        mutate(distance = round(distance, 4),
              diff = round(diff, 4),
              corr_dist = round(corr_dist, 4)) %>%
        # Arrange by shortest distance. If same, arrange by maximum phenotypic difference
        arrange(distance, desc(diff))

      # Initialize results variables
      selected_pairs <- data.frame()

      # Find the top pair based on distance
      top_pair <- distance_df %>%
        slice_head(n = 1) %>%
        # Create cluster column and give it a value of 1
        mutate(cluster = 1)

      # Add top pair to the reference dataframe
      selected_pairs <- rbind(selected_pairs, top_pair)

      # Update the list of selected species
      selected_species <- data.frame(top = top_pair$species1, bottom = top_pair$species2, diff = top_pair$diff, cluster = top_pair$cluster)

      # Establish the cluster matrix for calculating the Dunn index
      mat <- as.matrix(distance_matrix)

      # Initialize results for Dunn index
      dunn_results <- data.frame(species1 = character(), species2 = character(), Dunn_index = numeric(), diff = numeric(), cluster = numeric())
      dunn_result_cummulative <- data.frame()

      # Initialize variables
      current_dunn_index <- Inf  # Start with a very high value to ensure the loop runs
      iteration_limit <- 100
      iteration_count <- 0
      
      # Start loop to find the most convergent-divergent pairs
      while (current_dunn_index >= 1 && iteration_count < iteration_limit) {
        iteration_count <- iteration_count + 1
        query_cluster <- iteration_count + 1 # Twice the iteration count

        # Evaluate Dunn index for all candidate pairs
        dunn_candidates <- apply(distance_df, 1, function(row) {
          # Initialize cluster matrix
          query_species1 <- row["species1"]
          query_species2 <- row["species2"]
          query_diff <- row["diff"]

          # Dynamically assign selected species clusters
          selected_species <- rbind(selected_species, data.frame(top = query_species1, bottom = query_species2, diff = query_diff, cluster = query_cluster))

          # Build the cluster as a vector of integers, using as objects the cluster definitions and as names the species
          selected_species.str <- selected_species %>%
            # Pivot species in top and bottom to a single column
            pivot_longer(cols = c(top, bottom),
                        names_to = "label",
                        values_to = "species") %>%
            # Pull cluster number
            pull(cluster) %>%
            # Set names to species
            setNames(., selected_species %>%
                      pivot_longer(cols = c(top, bottom),
                                    names_to = "label",
                                    values_to = "species") %>%
                      pull(species))

          clusters <- as.integer(selected_species.str)

          # print(selected_species.str)
          # print(names(selected_species.str))

          # Rebuild distance matrix
          matrix <- mat[names(selected_species.str), names(selected_species.str)]
          dist_mat <- as.dist(matrix)

          # Compute Dunn index.
          # Use the modified version which computes intra-cluster distances only for a subset of clusters
          dunn_value <- mod_dunn(dist_mat, clusters, selected_cluster = query_cluster, verbose = FALSE)

          # Round Dunn index
          dunn_value <- round(dunn_value, 4)

          return(data.frame(species1 = query_species1, species2 = query_species2, Dunn_index = dunn_value, diff = query_diff, cluster = query_cluster))
        })

        # Combine results into a data frame
        dunn_candidates <- do.call(rbind, dunn_candidates)

        # Accumulate results for debug
        dunn_result_cummulative <- rbind(dunn_result_cummulative, dunn_candidates)
        dunn_result_cummulative <- dunn_result_cummulative %>%
          arrange(desc(cluster), desc(Dunn_index), desc(diff))

        # Select the top result with the highest Dunn index
        dunn_best <- dunn_candidates %>%
          arrange(desc(Dunn_index), desc(diff)) %>%
          slice_head(n = 1)

        # Update current Dunn index
        current_dunn_index <- dunn_best$Dunn_index

        # If the Dunn index is below 1, break the loop
        if (current_dunn_index < 1) {
          break
        }

        # Add the top result to cumulative results
        dunn_results <- rbind(dunn_results, dunn_best)

        # Update selected species
        selected_species <- rbind(selected_species, data.frame(top = dunn_best$species1, bottom = dunn_best$species2, diff = dunn_best$diff, cluster = dunn_best$cluster))
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

  outline = c(cycle.tag, fg.species.tag, bg.species.tag)
  simulated.traits.df[nrow(simulated.traits.df) +1,] <- outline
  
  chunk.counter = chunk.counter + 1
  
  # Write chunk to file when chunk size is reached or on final cycle
  if (chunk.counter > chunk.size || counter == as.integer(number.of.cycles)) {
    # Generate filename with zero-padded file number
    filename <- sprintf("resample_%03d.tab", file.counter)
    filepath <- file.path(outdir, filename)
    
    # Write the chunk
    write.table(simulated.traits.df, sep="\t", col.names=FALSE, row.names=FALSE, file=filepath, quote=FALSE)
    
    # Calculate progress and timing
    elapsed <- as.numeric(difftime(Sys.time(), start.time, units="secs"))
    pct.complete <- (counter / as.integer(number.of.cycles)) * 100
    cycles.per.sec <- counter / elapsed
    remaining.cycles <- as.integer(number.of.cycles) - counter
    eta.secs <- remaining.cycles / cycles.per.sec
    eta.mins <- eta.secs / 60
    
    logline <- sprintf("[%s] File %d: %s | Cycles %d-%d | Progress: %.1f%% | Elapsed: %.1f min | ETA: %.1f min",
                      format(Sys.time(), "%H:%M:%S"),
                      file.counter,
                      filename,
                      counter - nrow(simulated.traits.df) + 1,
                      counter,
                      pct.complete,
                      elapsed / 60,
                      eta.mins)
    write(logline, stdout())
    
    # Reset for next chunk
    simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))
    chunk.counter = 1
    file.counter = file.counter + 1
  }
}

# Final summary
end.time <- Sys.time()
total.elapsed <- as.numeric(difftime(end.time, start.time, units="mins"))
write(paste("[COMPLETE]", end.time, "|", number.of.cycles, "cycles in", round(total.elapsed, 2), "minutes"), stdout())
write(paste("[OUTPUT] Generated", file.counter - 1, "files in:", outdir), stdout())

