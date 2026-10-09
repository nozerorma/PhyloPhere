# phylo.R — Match the trait table to the species tree (pruning, clade MRCAs).
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Requires `trait_df`, `tree`, `tree_species`, `trait`, `resultsDir` (commons.R).
# The trait table and the tree come from NAME_CURATION, which settles the species names
# before this stage, so trait species are matched to tree tips by exact name (spaces read as
# underscores). Defines:
#   pruned_tree     the tree reduced to the species of trait_df
#   trait_df        reduced to the species of the tree; rows with an NA phenotype or
#                   count removed. Every row removed is reported with message()
#   find_taxon_mrca()  MRCA node and species count per taxon, for clade labels
# When 0.Data-pruning/ holds a pruned tree and trait table (written by
# 0.Data_pruning.Rmd), these are loaded instead of recomputed.
# =============================================================================

library(ape)

# Fallback logger, used only when commons.R has not defined debug_log().
if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# ── Tree and trait pruning ────────────────────────────────────────────────────

# Reuse the tree and trait table written by 0.Data_pruning.Rmd when present.
pruned_tree_path <- file.path(resultsDir, "0.Data-pruning", "pruned_tree_file.nwk")
pruned_trait_path <- file.path(resultsDir, "0.Data-pruning", "pruned_trait_file.tsv")

if (file.exists(pruned_tree_path) && file.exists(pruned_trait_path)) {
  debug_log("Found existing pruned tree and trait files in 0.Data-pruning. Loading them...")
  trait_df_ori <- trait_df # Table before pruning
  
  if (endsWith(pruned_trait_path, ".csv")) {
    p_sep <- ","
  } else if (endsWith(pruned_trait_path, ".tsv")) {
    p_sep <- "\t"
  } else {
    p_sep <- "\t"
  }
  
  trait_df <- read.csv(pruned_trait_path, sep = p_sep, stringsAsFactors = FALSE)
  trait_df[] <- lapply(trait_df, function(col) {
    if (is.character(col)) trimws(col) else col
  })
  
  pruned_tree <- ape::read.tree(file = pruned_tree_path)
  debug_log("Loaded pruned tree tips = %d, pruned trait rows = %d", length(pruned_tree$tip.label), nrow(trait_df))
} else {
  trait_df_ori <- trait_df # Table before pruning
  sp_norm <- gsub(" ", "_", trait_df$species)
  target <- ifelse(sp_norm %in% tree_species, sp_norm, NA_character_)
  discards <- sprintf("%s: no tree tip with this name", sp_norm[is.na(target)])

  # One row per tip: the first row in table order wins.
  rank_in_tip <- ave(seq_along(target), target, FUN = seq_along)
  dup <- !is.na(target) & rank_in_tip > 1
  discards <- c(discards, sprintf("%s: duplicate of the row kept for %s", sp_norm[dup], target[dup]))

  keep <- !is.na(target) & !dup
  trait_df$species <- target
  trait_df <- trait_df[keep, , drop = FALSE]
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, trait_df$species))
  debug_log("trait rows kept = %d, pruned tree tips = %d", nrow(trait_df), length(pruned_tree$tip.label))

  if (length(discards) > 0) {
    message(sprintf("[phylo.R] WARNING: %d trait row(s) removed while matching to the tree:\n  %s",
                    length(discards), paste(discards, collapse = "\n  ")))
  }
}

# ── Missing phenotypes ────────────────────────────────────────────────────────

# Drop species with NA in the primary trait (or in n_trait / c_trait when count
# data is in use), from the table and from the pruned tree.
if (exists("trait") && trait %in% names(trait_df)) {
  na_mask <- is.na(trait_df[[trait]])
  if (isTRUE(get0("has.n", ifnotfound = FALSE)) && exists("n_trait") && n_trait %in% names(trait_df)) {
    na_mask <- na_mask | is.na(trait_df[[n_trait]])
  }
  if (isTRUE(get0("has.c", ifnotfound = FALSE)) && exists("c_trait") && c_trait %in% names(trait_df)) {
    na_mask <- na_mask | is.na(trait_df[[c_trait]])
  }

  na_sp <- trait_df$species[na_mask]
  na_sp <- unique(na_sp[!is.na(na_sp)])

  if (length(na_sp) > 0) {
    debug_log("phylo.R: Removing %d species with NA phenotype/count values from dataset", length(na_sp))
    trait_df <- trait_df[!na_mask, , drop = FALSE]
    if (exists("pruned_tree") && !is.null(pruned_tree)) {
      pruned_tree <- ape::drop.tip(pruned_tree, intersect(pruned_tree$tip.label, na_sp))
    }
  }
}

# ── Clade MRCAs ───────────────────────────────────────────────────────────────

# Most recent common ancestor of each taxon, for the clade labels of the fan plots.
# `df` has one row per tree node with columns `taxa` (taxon of the node, NA when
# none), `node` and, optionally, `species`. Returns data.frame(taxa, mrca_node,
# n_species), with mrca_node NA when the taxon has no common ancestor among
# its nodes.
find_taxon_mrca <- function(tree, df) {
  library(ape)
  library(dplyr)
  
  # Path of node ids from `node` up to the root.
  trace_to_root <- function(node, tree) {
    path <- node
    current <- node
    edges <- tree$edge
    
    # Iteration cap, a guard against cycles in a malformed edge table.
    max_iter <- nrow(edges) + 10
    iter <- 0
    
    while(iter < max_iter) {
      iter <- iter + 1
      parent_row <- which(edges[, 2] == current)
      
      if(length(parent_row) == 0) {
        break
      }
      
      parent <- edges[parent_row[1], 1]  # Take first match if multiple
      path <- c(path, parent)
      current <- parent
    }
    
    return(path)
  }
  
  tip_labels <- tree$tip.label
  all_nodes <- c(1:length(tip_labels), unique(tree$edge[,1]))
  
  unique_taxa <- unique(df$taxa)
  unique_taxa <- unique_taxa[!is.na(unique_taxa)]
  
  result_taxa <- character()
  result_mrca <- integer()
  result_n_species <- integer()
  
  for(tax in unique_taxa) {
    taxon_data <- df[df$taxa == tax & !is.na(df$taxa), ]
    taxon_nodes <- taxon_data$node
    n_sp <- length(taxon_nodes)
    
    if(n_sp == 1) {
      mrca <- taxon_nodes[1]
    } else if(n_sp > 1) {
      # With species names matching tree tips, ape::getMRCA does the work.
      if("species" %in% colnames(df)) {
        taxon_species <- taxon_data$species
        
        tip_indices <- match(taxon_species, tree$tip.label)
        
        if(all(!is.na(tip_indices)) && length(tip_indices) > 1) {
          mrca <- getMRCA(tree, tip_indices)
        } else {
          # Otherwise the deepest node shared by the paths from each node to the root.
          all_paths <- lapply(taxon_nodes, trace_to_root, tree = tree)
          common_nodes <- Reduce(intersect, all_paths)
          
          if(length(common_nodes) > 0) {
            first_appearance <- sapply(common_nodes, function(cn) {
              min(which(all_paths[[1]] == cn))
            })
            mrca <- common_nodes[which.min(first_appearance)]
          } else {
            mrca <- NA_integer_
          }
        }
      } else {
        # No species column: same path-tracing rule.
        all_paths <- lapply(taxon_nodes, trace_to_root, tree = tree)
        common_nodes <- Reduce(intersect, all_paths)
        
        if(length(common_nodes) > 0) {
          first_appearance <- sapply(common_nodes, function(cn) {
            min(which(all_paths[[1]] == cn))
          })
          mrca <- common_nodes[which.min(first_appearance)]
        } else {
          mrca <- NA_integer_
        }
      }
    } else {
      mrca <- NA_integer_
    }
    
    result_taxa <- c(result_taxa, tax)
    result_mrca <- c(result_mrca, mrca)
    result_n_species <- c(result_n_species, n_sp)
  }
  
  result <- data.frame(
    taxa = result_taxa,
    mrca_node = result_mrca,
    n_species = result_n_species,
    stringsAsFactors = FALSE
  )
  
  return(result)
}
