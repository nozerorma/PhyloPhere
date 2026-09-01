# ----------------------------------------
# Phylogenetic Selection Algorithm for Independent Contrast Pairs
# ----------------------------------------
# Unified Phylogeny-Aware Pair Selection Engine
# 1. Calculates patristic distances between species
# 2. Selects species pairs that maximize trait differences while minimizing phylogenetic distance
# 3. Uses modified Dunn index to ensure phylogenetic independence between selected pairs
# ----------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Source the lean contrast selector engine
selector_script_path <- {
  this_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  this_dir <- if (!is.null(this_ofile)) dirname(this_ofile) else ""
  cand_paths <- c(
    file.path(getwd(), "src", "lean_contrast_selector.R"),
    file.path(getwd(), "..", "..", "CT", "local", "scripts", "lean_contrast_selector.R"),
    file.path(getwd(), "subworkflows", "CT", "local", "scripts", "lean_contrast_selector.R"),
    file.path(this_dir, "lean_contrast_selector.R"),
    file.path(this_dir, "..", "..", "CT", "local", "scripts", "lean_contrast_selector.R"),
    "/home/miguel/IBE-UPF/PhD/PhyloPhere/subworkflows/CT/local/scripts/lean_contrast_selector.R"
  )
  found <- cand_paths[nzchar(cand_paths) & file.exists(cand_paths)]
  if (length(found)) found[1] else ""
}

if (nzchar(selector_script_path) && file.exists(selector_script_path)) {
  source(selector_script_path)
} else {
  # Fallback inline definitions if path cannot be resolved
  mod_dunn_lean <- function(D, members, k) {
    c1 <- members[[k]]
    intra <- if (length(c1) > 1) max(D[c1, c1]) else 0
    if (intra == 0) return(Inf)
    inter <- Inf
    for (j in seq_along(members)) {
      if (j == k) next
      inter <- min(inter, min(D[c1, members[[j]]]))
    }
    if (!is.finite(inter)) return(Inf)
    inter / intra
  }
  overall_dunn_lean <- function(D, members) {
    if (length(members) <= 1) return(Inf)
    min(vapply(seq_along(members), function(k) mod_dunn_lean(D, members, k), numeric(1)))
  }
}

# Calculate patristic distances between species in a phylogenetic tree
calculate_patristic_distances <- function(tree, df) {
  species_list <- df$species
  debug_log("calculate_patristic_distances species = %d", length(species_list))
  
  species_in_tree <- intersect(species_list, tree$tip.label)
  if (length(species_in_tree) == 0) {
    stop("No species from data found in phylogenetic tree")
  }
  
  if (length(species_in_tree) < length(species_list)) {
    missing <- setdiff(species_list, tree$tip.label)
    warning(paste("Species not in tree:", paste(missing, collapse = ", ")))
    debug_log("calculate_patristic_distances missing in tree = %d", length(missing))
  }
  
  subtree <- ape::drop.tip(tree, setdiff(tree$tip.label, species_in_tree))
  debug_log("calculate_patristic_distances subtree tips = %d", length(subtree$tip.label))
  
  as.data.frame(as.matrix(ape::cophenetic.phylo(subtree)))
}

# Backward-compatible mod_dunn wrapper
mod_dunn <- function(distance = NULL, clusters, Data = NULL, method = "euclidean", 
                     selected_cluster, verbose = FALSE) {
  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method = method))
  if ("dist" %in% class(distance)) distance <- as.matrix(distance)
  
  sp_names <- rownames(distance)
  members <- lapply(sort(unique(clusters)), function(cl) sp_names[clusters == cl])
  mod_dunn_lean(distance, members, selected_cluster)
}

# Unified Canonical Contrast Selection
pair_sel.f <- function(distance_matrix, overlap_df, traits_df, my_trait) {
  mat <- as.matrix(distance_matrix)
  
  trait_df <- traits_df %>% dplyr::filter(trait == my_trait)
  if (has.n) {
    trait_df <- trait_df %>% dplyr::select(species, trait, value, n_data = n_trait)
  } else {
    trait_df <- trait_df %>% dplyr::select(species, trait, value)
  }

  distance_df <- as.data.frame(mat) %>%
    rownames_to_column(var = "species1") %>%
    gather(key = "species2", value = "distance", -species1) %>%
    dplyr::filter(distance != 0) %>%
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

  distance_df <- distance_df %>%
    left_join(overlap_df, by = c("species1", "species2")) %>%
    mutate(distance = round(distance, 4),
           diff = round(trait_diff, 4),
           abs_diff = round(abs(trait_diff), 4)) %>%
    dplyr::filter(!is.na(diff) & diff > 0)

  if (!"pss_score" %in% names(distance_df)) {
    distance_df$pss_score <- rep(NA_real_, nrow(distance_df))
  } else {
    distance_df$pss_score <- round(distance_df$pss_score, 4)
  }

  distance_df <- distance_df %>%
    {
      if (has.n) arrange(., distance, desc(abs_diff), desc(pair_n)) else arrange(., distance, desc(abs_diff))
    } %>%
    filter(!is.na(abs_diff)) %>%
    {
      cols <- c("species1", "species2", "distance", "abs_diff", "pss_score")
      if (has.n) cols <- c(cols, "pair_n")
      select(., any_of(cols))
    }

  debug_log("pair_sel.f candidate pairs = %d", nrow(distance_df))

  if (nrow(distance_df) == 0) {
    warning("pair_sel.f: no candidate contrast pairs. Returning 0 pairs.")
    empty_pairs <- data.frame(species1 = character(), species2 = character(), stringsAsFactors = FALSE)
    return(list(
      dunn_results = empty_pairs,
      selected_pairs = empty_pairs,
      dunn_result_cummulative = data.frame(),
      distance_df = distance_df,
      distance_matrix = distance_matrix
    ))
  }

  # Greedy selection with Dunn >= 1 gating
  selected_pairs <- data.frame()
  top_pair <- distance_df %>% slice_head(n = 1) %>% mutate(Dunn_index = Inf, cluster = 1)
  selected_pairs <- rbind(selected_pairs, top_pair)
  
  used_sp <- c(top_pair$species1, top_pair$species2)
  members <- list(used_sp)
  
  max_cap <- if (exists("max_contrasts", inherits = TRUE) &&
                !is.null(max_contrasts) &&
                !is.na(suppressWarnings(as.integer(max_contrasts))) &&
                as.integer(max_contrasts) > 0L) {
    as.integer(max_contrasts)
  } else if (exists("contrast_max_iter", inherits = TRUE) &&
             !is.null(contrast_max_iter) &&
             !is.na(suppressWarnings(as.integer(contrast_max_iter))) &&
             as.integer(contrast_max_iter) > 0L) {
    as.integer(contrast_max_iter) + 1L
  } else {
    Inf
  }

  dunn_results <- data.frame(
    species1 = character(), species2 = character(),
    Dunn_index = numeric(), abs_diff = numeric(),
    cluster = numeric(), stringsAsFactors = FALSE
  )
  if (has.n) dunn_results$pair_n <- numeric()

  repeat {
    if (length(members) >= max_cap) break
    
    avail <- !(distance_df$species1 %in% used_sp | distance_df$species2 %in% used_sp)
    if (!any(avail)) break
    
    cand <- distance_df[avail, , drop = FALSE]
    
    # Calculate modified Dunn index for each available candidate
    dunn_scores <- vapply(seq_len(nrow(cand)), function(i) {
      test_members <- c(members, list(c(cand$species1[i], cand$species2[i])))
      mod_dunn_lean(mat, test_members, length(test_members))
    }, numeric(1))
    
    cand$Dunn_index <- round(dunn_scores, 4)
    cand_passing <- cand[cand$Dunn_index >= 1.0, , drop = FALSE]
    if (nrow(cand_passing) == 0) break
    
    # Pick candidate that maximizes Dunn index (with trait diff and pair_n tie-breakers)
    best_cand <- cand_passing %>%
      {
        if (has.n) arrange(., desc(Dunn_index), desc(abs_diff), desc(pair_n))
        else arrange(., desc(Dunn_index), desc(abs_diff))
      } %>%
      slice_head(n = 1)
    
    next_cluster <- length(members) + 1L
    best_cand$cluster <- next_cluster
    
    # Verify overall Dunn across all clusters remains >= 1.0
    new_members <- c(members, list(c(best_cand$species1, best_cand$species2)))
    if (overall_dunn_lean(mat, new_members) < 1.0) break
    
    members <- new_members
    used_sp <- c(used_sp, best_cand$species1, best_cand$species2)
    selected_pairs <- rbind(selected_pairs, best_cand)
    dunn_results <- rbind(dunn_results, best_cand)
  }

  list(
    dunn_results = dunn_results,
    selected_pairs = selected_pairs,
    dunn_result_cummulative = dunn_results,
    distance_df = distance_df,
    distance_matrix = distance_matrix
  )
}

# FOP Parallel Multi-Hypothesis Selection with Voronoi Domain Partitioning
fop_pair_sel.f <- function(distance_matrix, overlap_df, traits_df, my_trait, max_fop = 100, seed = 42) {
  set.seed(seed)
  mat <- as.matrix(distance_matrix)
  
  # Step 1: Canonical greedy baseline (H1)
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
  
  # Step 4: Harvest non-replacement parallel hypotheses
  hypotheses <- list(H1 = canon_pairs)
  seen_sigs <- list(paste(sort(c(canon_pairs$species1, canon_pairs$species2)), collapse = "|"))
  summary_rows <- list()
  
  canon_dunn <- overall_dunn_lean(mat, canon_members)
  h1_mean_pss <- if ("pss_score" %in% names(canon_pairs)) round(mean(canon_pairs$pss_score, na.rm = TRUE), 4) else NA_real_
  summary_rows[[1]] <- data.frame(
    hypothesis_id = "H1",
    is_canonical = TRUE,
    num_pairs = K,
    min_dunn = round(canon_dunn, 4),
    mean_distance = round(mean(canon_pairs$distance), 4),
    mean_abs_diff = round(mean(canon_pairs$abs_diff), 4),
    mean_pss_score = h1_mean_pss,
    jaccard_to_h1 = 1.0,
    pair_composition = paste(paste(canon_pairs$species1, canon_pairs$species2, sep = "~"), collapse = "; "),
    stringsAsFactors = FALSE
  )
  
  # Combinatorial harvest across Voronoi domains
  pool_sizes <- sapply(alt_pools, nrow)
  total_combos <- prod(pmax(pool_sizes, 1))
  
  if (total_combos > 1) {
    n_sample <- min(as.numeric(total_combos), max_fop * 20)
    
    for (iter in seq_len(n_sample)) {
      if (length(hypotheses) >= max_fop) break
      
      chosen_idx <- sapply(seq_len(K), function(k) {
        if (pool_sizes[k] == 0) return(NA_integer_)
        sample.int(pool_sizes[k], 1)
      })
      
      if (any(is.na(chosen_idx))) next
      
      cand_pairs_list <- lapply(seq_len(K), function(k) alt_pools[[k]][chosen_idx[k], ])
      cand_h_df <- do.call(rbind, cand_pairs_list)
      
      all_sp <- c(cand_h_df$species1, cand_h_df$species2)
      if (length(unique(all_sp)) < 2 * K) next
      
      sig <- paste(sort(all_sp), collapse = "|")
      if (sig %in% seen_sigs) next
      
      h_members <- lapply(seq_len(K), function(i) c(cand_h_df$species1[i], cand_h_df$species2[i]))
      h_dunn <- overall_dunn_lean(mat, h_members)
      
      if (h_dunn >= 1.0) {
        h_id <- paste0("H", length(hypotheses) + 1)
        cand_h_df$cluster <- seq_len(K)
        hypotheses[[h_id]] <- cand_h_df
        seen_sigs <- c(seen_sigs, sig)
        
        h1_sp <- c(canon_pairs$species1, canon_pairs$species2)
        jaccard <- length(intersect(h1_sp, all_sp)) / length(union(h1_sp, all_sp))
        h_mean_pss <- if ("pss_score" %in% names(cand_h_df)) round(mean(cand_h_df$pss_score, na.rm = TRUE), 4) else NA_real_
        
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          hypothesis_id = h_id,
          is_canonical = FALSE,
          num_pairs = K,
          min_dunn = round(h_dunn, 4),
          mean_distance = round(mean(cand_h_df$distance), 4),
          mean_abs_diff = round(mean(cand_h_df$abs_diff), 4),
          mean_pss_score = h_mean_pss,
          jaccard_to_h1 = round(jaccard, 4),
          pair_composition = paste(paste(cand_h_df$species1, cand_h_df$species2, sep = "~"), collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  summary_df <- do.call(rbind, summary_rows)
  
  list(
    canon_pairs = canon_pairs,
    hypotheses = hypotheses,
    summary_df = summary_df
  )
}
