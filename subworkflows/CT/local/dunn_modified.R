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
