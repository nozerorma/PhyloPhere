# ----------------------------------------
# Phylogenetic Selection Algorithm for Independent Contrast Pairs
# ----------------------------------------
# 
# This module implements a phylogeny-aware pair selection algorithm that:
# 1. Calculates patristic (evolutionary) distances between species
# 2. Selects species pairs that maximize trait differences while minimizing phylogenetic distance
# 3. Uses a modified Dunn index to ensure phylogenetic independence between selected pairs
#
# The algorithm ensures that selected pairs are:
# - Phylogenetically close within each pair (small patristic distance)
# - Phenotypically divergent within each pair (large trait difference)
# - Phylogenetically independent between pairs (Dunn index ≥ 1)
# ----------------------------------------
# This is my most original contribution to PhyloPhere! Heheh

# Load required libraries
library(ape)
library(dplyr)
library(tidyr)

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Calculate patristic distances between species in a phylogenetic tree
#
# Args:
#   tree: phylo object (full or pre-pruned phylogenetic tree)
#   df: Data frame containing a 'species' column with species names
#
# Returns:
#   Data frame (matrix format) with pairwise patristic distances


calculate_patristic_distances <- function(tree, df) {
  species_list <- df$species
  debug_log("calculate_patristic_distances species = %d", length(species_list))
  
  # Prune tree to include only species present in the data
  species_in_tree <- intersect(species_list, tree$tip.label)
  
  if (length(species_in_tree) == 0) {
    stop("No species from data found in phylogenetic tree")
  }
  
  if (length(species_in_tree) < length(species_list)) {
    missing <- setdiff(species_list, tree$tip.label)
    warning(paste("Species not in tree:", paste(missing, collapse = ", ")))
    debug_log("calculate_patristic_distances missing in tree = %d", length(missing))
  }
  
  # Prune tree to selected species
  subtree <- ape::drop.tip(tree, setdiff(tree$tip.label, species_in_tree))
  debug_log("calculate_patristic_distances subtree tips = %d", length(subtree$tip.label))
  
  # Calculate patristic distances
  patristic_distances <- ape::cophenetic.phylo(subtree)
  
  # Convert to data frame format
  distances_df <- as.data.frame(as.matrix(patristic_distances))
  
  return(distances_df)
}

# Modified Dunn index for phylogenetic cluster validation
#
# Calculates the Dunn index for a specific cluster to assess separation from other clusters.
# Used to ensure phylogenetic independence between selected species pairs.
#
# Dunn index = min(inter-cluster distance) / max(intra-cluster distance)
# Values ≥ 1 indicate well-separated clusters (suitable for independent contrasts)
#
# Args:
#   distance: Distance matrix (or dist object) between all points
#   clusters: Integer vector assigning each point to a cluster
#   Data: Optional data matrix to compute distances from (if distance not provided)
#   method: Distance metric for computing distances from Data
#   selected_cluster: ID of the cluster to evaluate
#   verbose: If TRUE, print diagnostic messages
#
# Returns:
#   Numeric Dunn index value for the selected cluster

mod_dunn <- function(distance = NULL, clusters, Data = NULL, method = "euclidean", 
                     selected_cluster, verbose = FALSE) {
  
  # Ensure either 'distance' or 'Data' is provided
  if (is.null(distance) & is.null(Data)) {
    stop("One of 'distance' or 'Data' is required")
  }
  
  # Compute distance matrix if not provided
  if (is.null(distance)) {
    distance <- as.matrix(dist(Data, method = method))
    if (verbose) message("Distance matrix computed from Data using 'dist' function.")
  }
  
  # Convert 'dist' object to matrix if necessary
  if ("dist" %in% class(distance)) {
    distance <- as.matrix(distance)
    if (verbose) message("Distance object converted to a matrix.")
  }
  
  # Validate 'selected_cluster'
  nc <- max(clusters)
  if (!(selected_cluster %in% seq_len(nc))) {
    stop("Invalid cluster ID in 'selected_cluster'")
  }
  
  # Initialize variables for distances
  interClust <- rep(NA, nc) # Inter-cluster distances for selected cluster
  intraClust <- NA          # Intra-cluster distance for selected cluster
  
  # Get indices for the selected cluster
  c1 <- which(clusters == selected_cluster)
  
  # Calculate intra-cluster distance for the selected cluster
  if (length(c1) > 1) {
    intraClust <- max(distance[c1, c1])
  } else {
    intraClust <- 0 # Single-point cluster
  }
  if (verbose) message(paste("Intra-cluster distance for cluster", selected_cluster, "calculated:", intraClust))
  
  # Calculate inter-cluster distances only between the selected cluster and the rest
  for (j in seq_len(nc)) {
    if (j == selected_cluster) next # Skip self-comparison
    
    c2 <- which(clusters == j) # Indices of points in cluster j
    interClust[j] <- min(distance[c1, c2])
    if (verbose) message(paste("Inter-cluster distance between cluster", selected_cluster, "and cluster", j, "calculated:", interClust[j]))
  }
  
  # Calculate Dunn index
  interDist <- min(interClust, na.rm = TRUE)
  
  if (intraClust == 0) {
    if (verbose) message("Warning: Intra-cluster distance is zero, Dunn index may not be meaningful.")
    dunn <- Inf
  } else {
    dunn <- interDist / intraClust
  }
  
  if (verbose) message(paste("Dunn index calculated:", dunn))
  
  return(dunn)
}

# Species pair selection algorithm for independent phylogenetic contrasts
#
# Selects species pairs that:
# 1. Have small phylogenetic distances (close evolutionary relationship)
# 2. Have large trait value differences (phenotypic divergence)
# 3. Maintain phylogenetic independence between pairs (Dunn index ≥ 1)
#
# The algorithm iteratively adds pairs, evaluating each candidate using the Dunn index
# to ensure new pairs don't overlap phylogenetically with previously selected pairs.
#
# Args:
#   distance_matrix: Patristic distance matrix from calculate_patristic_distances()
#   overlap_df: Data frame with columns 'species1', 'species2', 'absdiff' (trait difference)
#   traits_df: Data frame with columns 'species', 'trait', 'value'
#   my_trait: Name of the trait to use for pair selection
#
# Returns:
#   List containing:
#     - dunn_results: Pairs selected with Dunn ≥ 1
#     - selected_pairs: All selected pairs (first pair + Dunn-validated pairs)
#     - dunn_result_cummulative: All candidate pairs evaluated with Dunn indices
#     - distance_df: Processed distance dataframe with trait differences
#     - distance_matrix: Original distance matrix

pair_sel.f <- function(distance_matrix, overlap_df, traits_df, my_trait) {
  # Transform distance matrix into long-format dataframe with trait values
  # Keep optional sampling column `n` if provided to support tertiary tie-breaking
  if (has.p) {
    debug_log("pair_sel.f has population data")
  } else {
    debug_log("pair_sel.f does not have population data")
  }

  trait_df <- traits_df %>%
  
  if (has.p) {
    trait_df <- trait_df %>%
      dplyr::filter(trait == my_trait) %>%
      dplyr::select(species, trait, p_data = p_trait)
  } else {
    trait_df <- trait_df %>%
      dplyr::filter(trait == my_trait) %>%
      dplyr::select(species, trait)
  }


  distance_df <- as.data.frame(as.matrix(distance_matrix)) %>%
    rownames_to_column(var = "species1") %>%
    gather(key = "species2", value = "distance", -species1) %>%
    dplyr::filter(distance != 0) %>%  # Remove self-comparisons
    left_join(trait_df, by = c("species1" = "species")) %>%
    dplyr::rename(value1 = value) %>%
    left_join(trait_df, by = c("species2" = "species")) %>%
    dplyr::rename(value2 = value) %>%
    {
      if (has.p) {
        dplyr::rename(., p1 = p_data.x, p2 = p_data.y) %>%
          dplyr::mutate(pair_p = ifelse(is.na(p1) | is.na(p2), NA_real_, p1 + p2))
      } else {
        .
      }
    }

  # Filter to valid pairs from overlap analysis and format
  distance_df <- distance_df %>%
    left_join(overlap_df, by = c("species1", "species2")) %>%
    mutate(distance = round(distance, 4),
           diff = round(abs(trait_diff), 4)) %>%
    {
      if (has.p) arrange(., distance, desc(diff), desc(pair_p)) else arrange(., distance, desc(diff))
    } %>%
    filter(!is.na(diff)) %>%
    {
      if (has.p) select(., species1, species2, distance, diff, pair_p) else select(., species1, species2, distance, diff)
    }
  debug_log("pair_sel.f candidate pairs = %d", nrow(distance_df))

  # Select initial pair: smallest phylogenetic distance with largest trait difference
  selected_pairs <- data.frame()

  top_pair <- distance_df %>%
    {
      if (has.p) arrange(., distance, desc(diff), desc(pair_p)) else arrange(., distance, desc(diff))
    } %>%
    slice_head(n = 1) %>%
    mutate(cluster = 1)

  selected_pairs <- rbind(selected_pairs, top_pair)
  debug_log("pair_sel.f top pair = %s vs %s", top_pair$species1, top_pair$species2)
  selected_species <- data.frame(
    top = top_pair$species1, 
    bottom = top_pair$species2, 
    diff = top_pair$diff, 
    cluster = top_pair$cluster
  )
  mat <- as.matrix(distance_matrix)

  # Initialize tracking variables
  dunn_results <- data.frame(
    species1 = character(), 
    species2 = character(), 
    Dunn_index = numeric(), 
    diff = numeric(), 
    cluster = numeric(),
    stringsAsFactors = FALSE
  )
  if (has.p) dunn_results$pair_p <- numeric()
  dunn_result_cummulative <- data.frame()

  current_dunn_index <- Inf
  iteration_limit <- 10
  iteration_count <- 0

  # Iteratively add pairs while maintaining phylogenetic independence (Dunn ≥ 1)
  while (current_dunn_index >= 1 && iteration_count < iteration_limit) {
    iteration_count <- iteration_count + 1
    querry_cluster <- iteration_count + 1

    message(paste0("Pair selection iteration: ", iteration_count))

    # Evaluate all candidate pairs for Dunn index
    dunn_candidates <- apply(distance_df, 1, function(row) {
      query_species1 <- row["species1"]
      query_species2 <- row["species2"]
      query_diff <- row["diff"]
      query_pair_p <- if (has.p) row["pair_p"] else NA_real_

      # Temporarily add candidate pair as new cluster
      selected_species <- rbind(
        selected_species, 
        data.frame(
          top = query_species1, 
          bottom = query_species2, 
          diff = query_diff, 
          cluster = querry_cluster
        )
      )

      # Create cluster assignment vector
      selected_species.str <- selected_species %>%
        pivot_longer(cols = c(top, bottom), names_to = "label", values_to = "species") %>%
        pull(cluster) %>%
        setNames(., selected_species %>%
                   pivot_longer(cols = c(top, bottom), names_to = "label", values_to = "species") %>%
                   pull(species))

      clusters <- as.integer(selected_species.str)
      matrix <- mat[names(selected_species.str), names(selected_species.str)]
      dist_mat <- as.dist(matrix)

      # Calculate Dunn index for candidate cluster
      dunn_value <- mod_dunn(
        dist_mat, 
        clusters, 
        selected_cluster = querry_cluster, 
        verbose = FALSE
      )
      dunn_value <- round(dunn_value, 4)

      out <- data.frame(
        species1 = query_species1, 
        species2 = query_species2, 
        Dunn_index = dunn_value, 
        diff = query_diff, 
        cluster = querry_cluster,
        stringsAsFactors = FALSE
      )

      if (has.p) out$pair_p <- as.numeric(query_pair_p)

      return(out)
    })

    dunn_candidates <- do.call(rbind, dunn_candidates)
    dunn_result_cummulative <- rbind(dunn_result_cummulative, dunn_candidates) %>%
      {
        if (has.p) arrange(., desc(cluster), desc(Dunn_index), desc(diff), desc(pair_p)) else arrange(., desc(cluster), desc(Dunn_index), desc(diff))
      }

    # Select best candidate: highest Dunn index, then largest trait difference
    dunn_best <- dunn_candidates %>%
      {
        if (has.p) arrange(., desc(Dunn_index), desc(diff), desc(pair_p)) else arrange(., desc(Dunn_index), desc(diff))
      } %>%
      slice_head(n = 1)

    current_dunn_index <- dunn_best$Dunn_index

    # Stop if best candidate doesn't maintain phylogenetic independence
    if (current_dunn_index < 1) {
      message(paste0("Stopping: Best candidate Dunn index (", 
                     round(current_dunn_index, 4), 
                     ") < 1. Selected ", 
                     iteration_count, 
                     " phylogenetically independent pairs."))
      break
    }

    dunn_results <- rbind(dunn_results, dunn_best)
    selected_species <- rbind(
      selected_species, 
      data.frame(
        top = dunn_best$species1, 
        bottom = dunn_best$species2, 
        diff = dunn_best$diff, 
        cluster = dunn_best$cluster
      )
    )
  }

  # Combine initial pair with Dunn-validated pairs
  selected_pairs <- selected_pairs %>%
    dplyr::select(species1, species2)

  dunn_out <- dunn_results %>%
    dplyr::select(species1, species2)

  selected_pairs <- bind_rows(selected_pairs, dunn_out)

  if (iteration_count == iteration_limit) {
    warning(paste0("Iteration limit (", iteration_limit, 
                   ") reached before Dunn index dropped below 1. ",
                   "Consider increasing iteration_limit if more pairs are expected."))
  }

  return(list(
    dunn_results = dunn_results,
    selected_pairs = selected_pairs,
    dunn_result_cummulative = dunn_result_cummulative,
    distance_df = distance_df,
    distance_matrix = distance_matrix
  ))
}
