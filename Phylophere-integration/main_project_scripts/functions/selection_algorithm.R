# Define function to calculate patristic distances
calculate_patristic_distances <- function(tree, df) {
  species_list <- df$species  # Extract species list from dataframe
  
  # Subset tree to include only species in the list
  subtree <- drop.tip(tree, setdiff(tree$tip.label, species_list))

  # Calculate patristic distances between species
  patristic_distances <- cophenetic(subtree)
  
  # Convert distance matrix to dataframe
  distances_df <- as.data.frame(as.matrix(patristic_distances))
  colnames(distances_df) <- rownames(distances_df)  # Set column names to species
  
  return(distances_df)
}

# Modified Dunn index function
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

# Species selection algorithm
pair_sel.f <- function(distance_matrix, overlap_df, traits_df, my_trait) {
  # Transform the distance matrix into a dataframe
  trait_df <- traits_df %>%
    filter(trait == my_trait) %>%
    dplyr::select(species, value)

  distance_df <- as.data.frame(as.matrix(distance_matrix)) %>%
    rownames_to_column(var = "species1") %>%
    gather(key = "species2", value = "distance", -species1) %>%
    filter(distance != 0) %>%
    left_join(trait_df, by = c("species1" = "species")) %>%
    dplyr::rename(value1 = value) %>%
    left_join(trait_df, by = c("species2" = "species")) %>%
    dplyr::rename(value2 = value)

  # Filter out pairs not in the overlap dataframe
  distance_df <- distance_df %>%
    left_join(overlap_df, by = c("species1", "species2")) %>%
    mutate(distance = round(distance, 4),
           diff = round(abs(absdiff), 4)) %>%
    arrange(distance, desc(diff)) %>%
    filter(!is.na(diff)) %>%
    select(species1, species2, distance, diff)

  # Initialize results variables
  selected_pairs <- data.frame()

  # Find the top pair based on distance
  top_pair <- distance_df %>%
    arrange(distance) %>%
    slice_head(n = 1) %>%
    mutate(cluster = 1)

  selected_pairs <- rbind(selected_pairs, top_pair)
  selected_species <- data.frame(top = top_pair$species1, bottom = top_pair$species2, diff = top_pair$diff, cluster = top_pair$cluster)
  mat <- as.matrix(distance_matrix)

  # Initialize results for Dunn index
  dunn_results <- data.frame(species1 = character(), species2 = character(), Dunn_index = numeric(), diff = numeric(), cluster = numeric())
  dunn_result_cummulative <- data.frame()

  # Initialize variables
  current_dunn_index <- Inf
  iteration_limit <- 10
  iteration_count <- 0

  # Begin loop
  while (current_dunn_index >= 1 && iteration_count < iteration_limit) {
    iteration_count <- iteration_count + 1
    querry_cluster <- iteration_count + 1

    print(paste0("Iteration: ", iteration_count))

    dunn_candidates <- apply(distance_df, 1, function(row) {
      query_species1 <- row["species1"]
      query_species2 <- row["species2"]
      query_diff <- row["diff"]

      selected_species <- rbind(selected_species, data.frame(top = query_species1, bottom = query_species2, diff = query_diff, cluster = querry_cluster))

      selected_species.str <- selected_species %>%
        pivot_longer(cols = c(top, bottom), names_to = "label", values_to = "species") %>%
        pull(cluster) %>%
        setNames(., selected_species %>%
                   pivot_longer(cols = c(top, bottom), names_to = "label", values_to = "species") %>%
                   pull(species))

      clusters <- as.integer(selected_species.str)
      matrix <- mat[names(selected_species.str), names(selected_species.str)]
      dist_mat <- as.dist(matrix)

      dunn_value <- mod_dunn(dist_mat, clusters, selected_cluster = querry_cluster, verbose = FALSE)
      dunn_value <- round(dunn_value, 4)

      return(data.frame(species1 = query_species1, species2 = query_species2, Dunn_index = dunn_value, diff = query_diff, cluster = querry_cluster))
    })

    dunn_candidates <- do.call(rbind, dunn_candidates)
    dunn_result_cummulative <- rbind(dunn_result_cummulative, dunn_candidates) %>%
      arrange(desc(cluster), desc(Dunn_index), desc(diff))

    dunn_best <- dunn_candidates %>%
      arrange(desc(Dunn_index), desc(diff)) %>%
      slice_head(n = 1)

    current_dunn_index <- dunn_best$Dunn_index

    if (current_dunn_index < 1) {
      break
    }

    dunn_results <- rbind(dunn_results, dunn_best)
    selected_species <- rbind(selected_species, data.frame(top = dunn_best$species1, bottom = dunn_best$species2, diff = dunn_best$diff, cluster = dunn_best$cluster))
  }

  selected_pairs <- selected_pairs %>%
    dplyr::select(species1, species2)

  dunn_out <- dunn_results %>%
    dplyr::select(species1, species2)

  selected_pairs <- bind_rows(selected_pairs, dunn_out)

  if (iteration_count == iteration_limit) {
    warning("Iteration limit reached before Dunn index dropped below 1.")
  }

  return(list(dunn_results = dunn_results,
              selected_pairs = selected_pairs,
              dunn_result_cummulative = dunn_result_cummulative,
              distance_df = distance_df,
              distance_matrix = distance_matrix))
}
