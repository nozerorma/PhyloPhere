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
library(tibble)

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
  if (has.n) {
    debug_log("pair_sel.f has population data")
  } else {
    debug_log("pair_sel.f does not have population data")
  }

  trait_df <- traits_df
  if (has.n) {
    trait_df <- trait_df %>%
      dplyr::filter(trait == my_trait) %>%
      dplyr::select(species, trait, value, n_data = n_trait)
  } else {
    trait_df <- trait_df %>%
      dplyr::filter(trait == my_trait) %>%
      dplyr::select(species, trait, value)
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
      if (has.n) {
        dplyr::rename(., n1 = n_data.x, n2 = n_data.y) %>%
          dplyr::mutate(pair_n = ifelse(is.na(n1) | is.na(n2), NA_real_, n1 + n2))
      } else {
        .
      }
    }

  # Filter to valid pairs from overlap analysis and format
  distance_df <- distance_df %>%
    left_join(overlap_df, by = c("species1", "species2")) %>%
    mutate(distance = round(distance, 4),
          diff = round(trait_diff, 4),
          abs_diff = round(abs(trait_diff), 4)) %>%
    # Filter out pairs with under 0 difff
    dplyr::filter(!is.na(diff) & diff > 0) %>%  
    {
      if (has.n) arrange(., distance, desc(abs_diff), desc(pair_n)) else arrange(., distance, desc(abs_diff))
    } %>%
    filter(!is.na(abs_diff)) %>%
    {
      if (has.n) select(., species1, species2, distance, abs_diff, pair_n) else select(., species1, species2, distance, abs_diff)
    }
  debug_log("pair_sel.f candidate pairs = %d", nrow(distance_df))

  # Select initial pair: smallest phylogenetic distance with largest trait difference
  selected_pairs <- data.frame()

  top_pair <- distance_df %>%
    {
      if (has.n) arrange(., distance, desc(abs_diff), desc(pair_n)) else arrange(., distance, desc(abs_diff))
    } %>%
    slice_head(n = 1) %>%
    mutate(cluster = 1)

  selected_pairs <- rbind(selected_pairs, top_pair)
  debug_log("pair_sel.f top pair = %s vs %s", top_pair$species1, top_pair$species2)
  selected_species <- data.frame(
    top = top_pair$species1, 
    bottom = top_pair$species2, 
    abs_diff = top_pair$abs_diff, 
    cluster = top_pair$cluster
  )
  mat <- as.matrix(distance_matrix)

  # Initialize tracking variables
  dunn_results <- data.frame(
    species1 = character(), 
    species2 = character(), 
    Dunn_index = numeric(), 
    abs_diff = numeric(), 
    cluster = numeric(),
    stringsAsFactors = FALSE
  )
  if (has.n) dunn_results$pair_n <- numeric()
  dunn_result_cummulative <- data.frame()

  current_dunn_index <- Inf
  # Allow override via contrast_max_iter global (set from params in commons.R)
  iteration_limit <- if (exists("contrast_max_iter", inherits = TRUE)) as.integer(contrast_max_iter) else 3L
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
      query_abs_diff <- row["abs_diff"]
      query_pair_n <- if (has.n) row["pair_n"] else NA_real_

      # Temporarily add candidate pair as new cluster
      selected_species <- rbind(
        selected_species, 
        data.frame(
          top = query_species1, 
          bottom = query_species2, 
          abs_diff = query_abs_diff, 
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
        abs_diff = query_abs_diff, 
        cluster = querry_cluster,
        stringsAsFactors = FALSE
      )

      if (has.n) out$pair_n <- as.numeric(query_pair_n)

      return(out)
    })

    dunn_candidates <- do.call(rbind, dunn_candidates)
    dunn_result_cummulative <- rbind(dunn_result_cummulative, dunn_candidates) %>%
      {
        if (has.n) arrange(., desc(cluster), desc(Dunn_index), desc(abs_diff), desc(pair_n)) else arrange(., desc(cluster), desc(Dunn_index), desc(abs_diff))
      }

    # Select best candidate: highest Dunn index, then largest trait difference
    dunn_best <- dunn_candidates %>%
      {
        if (has.n) arrange(., desc(Dunn_index), desc(abs_diff), desc(pair_n)) else arrange(., desc(Dunn_index), desc(abs_diff))
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
        abs_diff = dunn_best$abs_diff, 
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

# Overall Dunn index across all clusters in a hypothesis set
overall_dunn_calc <- function(mat, members) {
  if (length(members) <= 1) return(Inf)
  dunn_vals <- numeric(length(members))
  for (k in seq_along(members)) {
    c1 <- members[[k]]
    intra <- if (length(c1) > 1) max(mat[c1, c1]) else 0
    if (intra == 0) {
      dunn_vals[k] <- Inf
      next
    }
    inter <- Inf
    for (j in seq_along(members)) {
      if (j == k) next
      c2 <- members[[j]]
      inter <- min(inter, min(mat[c1, c2]))
    }
    dunn_vals[k] <- if (is.finite(inter)) inter / intra else Inf
  }
  min(dunn_vals)
}

# FOP Non-Greedy Pair Selection Algorithm
# Generates canonical baseline (H1) and harvests parallel independent hypotheses (H2...HM)
fop_pair_sel.f <- function(distance_matrix, overlap_df, traits_df, my_trait, max_fop = 100, seed = 42) {
  set.seed(seed)
  mat <- as.matrix(distance_matrix)
  
  # Step 1: Run canonical greedy selection
  canon_res <- pair_sel.f(distance_matrix, overlap_df, traits_df, my_trait)
  canon_pairs <- canon_res$selected_pairs
  K <- nrow(canon_pairs)
  
  if (K == 0) {
    return(list(
      canon_pairs = canon_pairs,
      hypotheses = list(),
      summary_df = data.frame()
    ))
  }
  
  canon_members <- lapply(seq_len(K), function(i) c(canon_pairs$species1[i], canon_pairs$species2[i]))
  
  # Step 2: Voronoi Clade Domain Partitioning
  all_species <- rownames(mat)
  species_domain <- setNames(integer(length(all_species)), all_species)
  for (sp in all_species) {
    dists_to_canon <- sapply(seq_len(K), function(k) {
      min(mat[sp, canon_members[[k]][1]], mat[sp, canon_members[[k]][2]])
    })
    species_domain[sp] <- which.min(dists_to_canon)
  }
  
  # Step 3: Domain Pools from valid candidates
  cand_df <- canon_res$distance_df
  alt_pools <- list()
  for (k in seq_len(K)) {
    dom_sp <- names(species_domain)[species_domain == k]
    alt_pools[[k]] <- cand_df %>%
      dplyr::filter(species1 %in% dom_sp & species2 %in% dom_sp)
  }
  
  # Step 4: Harvest Hypotheses (H1 = Canonical, plus random draws without replacement)
  harvested_list <- list()
  unique_sigs <- character()

  # H1 Canonical Baseline
  h1_members <- canon_members
  h1_sig <- paste(sapply(h1_members, function(x) paste(sort(x), collapse = "-")), collapse = " | ")
  h1_dunn <- overall_dunn_calc(mat, h1_members)
  h1_sum_diff <- sum(sapply(seq_len(K), function(i) {
    sub <- cand_df %>% dplyr::filter(species1 == canon_pairs$species1[i] & species2 == canon_pairs$species2[i])
    if (nrow(sub) > 0) sub$abs_diff[1] else 0
  }))

  unique_sigs <- c(unique_sigs, h1_sig)
  harvested_list[[1]] <- list(
    members = h1_members,
    overall_dunn = h1_dunn,
    sum_diff = h1_sum_diff,
    id = "H1"
  )

  # Stochastic draws for alternative hypotheses
  n_draws <- 1000
  for (draw in seq_len(n_draws)) {
    if (length(harvested_list) >= max_fop) break
    shuffled_cands <- cand_df[sample(nrow(cand_df)), ]
    current_members <- list()
    used_sp <- character()

    for (r in seq_len(nrow(shuffled_cands))) {
      t <- shuffled_cands$top[r]; if (is.null(t)) t <- shuffled_cands$species1[r]
      b <- shuffled_cands$bot[r]; if (is.null(b)) b <- shuffled_cands$species2[r]
      if (t %in% used_sp || b %in% used_sp) next

      test_members <- c(current_members, list(c(t, b)))
      if (overall_dunn_calc(mat, test_members) >= 1.0) {
        current_members <- test_members
        used_sp <- c(used_sp, t, b)
      }
    }

    if (length(current_members) == K) {
      pair_strs <- sapply(current_members, function(x) paste(sort(x), collapse = "-"))
      sorted_sig <- paste(sort(pair_strs), collapse = " | ")

      if (!sorted_sig %in% unique_sigs) {
        unique_sigs <- c(unique_sigs, sorted_sig)
        ov_dunn <- overall_dunn_calc(mat, current_members)
        sum_diff <- sum(sapply(current_members, function(pair) {
          sub <- cand_df %>% dplyr::filter(species1 == pair[1] & species2 == pair[2])
          if (nrow(sub) > 0) sub$abs_diff[1] else 0
        }))

        harvested_list[[length(harvested_list) + 1]] <- list(
          members = current_members,
          overall_dunn = ov_dunn,
          sum_diff = sum_diff,
          id = paste0("H", length(harvested_list) + 1)
        )
      }
    }
  }

  # Build summary dataframe
  summary_rows <- list()
  hypotheses_dict <- list()

  for (idx in seq_along(harvested_list)) {
    h_item <- harvested_list[[idx]]
    h_id <- paste0("H", idx)

    # Format pair dataframe
    pair_df <- data.frame(
      species1 = sapply(h_item$members, `[`, 1),
      species2 = sapply(h_item$members, `[`, 2),
      stringsAsFactors = FALSE
    )
    hypotheses_dict[[h_id]] <- pair_df

    summary_rows[[idx]] <- data.frame(
      hypothesis_id = h_id,
      n_pairs = length(h_item$members),
      sum_abs_diff = round(h_item$sum_diff, 4),
      overall_dunn = round(h_item$overall_dunn, 4),
      is_canonical = (idx == 1),
      pairs_str = paste(sapply(h_item$members, function(p) paste0("(", p[1], " vs ", p[2], ")")), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)

  return(list(
    canon_pairs = canon_pairs,
    hypotheses = hypotheses_dict,
    summary_df = summary_df,
    domain_pools = alt_pools
  ))
}

