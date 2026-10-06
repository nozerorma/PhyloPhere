# phylo.R — Match the trait table to the species tree (tax_id mapping, pruning, clade MRCAs).
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Requires `trait_df`, `tree`, `tree_species`, `trait`, `resultsDir` (commons.R)
# and, optionally, `params$tax_id` (a CSV/TSV with `tax_id` and `species`
# columns, giving tree-side species names). Defines:
#   has.TAX_ID      TRUE when trait_df carries a usable `tax_id` column
#   pruned_tree     the tree reduced to the species of trait_df
#   trait_df        reduced to the species of the tree (tip names, via tax_id
#                   when available); rows with an NA phenotype or count removed
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

# ── tax_id mapping ────────────────────────────────────────────────────────────

# TAX_ID mode needs either a tax_id column in trait_df or a mapping file.
has.TAX_ID <- FALSE
tax_id_file <- if (exists("params")) params$tax_id else ""
debug_log("tax_id_file = %s", ifelse(nzchar(tax_id_file), tax_id_file, "<none>"))
# Separator from the file extension (comma when it is neither .csv nor .tsv).
if (endsWith(tax_id_file, ".csv")) {
  sep_char <- ","
} else if (endsWith(tax_id_file, ".tsv")) {
  sep_char <- "\t"
} else {
  sep_char <- ","
}

# Give a distinct tax_id to every species that shares one (e.g. several GenBank
# accessions of the same organism as distinct tree tips). Mirrors the synthetic
# tax_id assignment of CT_DISAMBIGUATION/local/src/phylo/species_mapping.py: the
# first species (alphabetically) keeps the original tax_id and every other one
# gets tax_id + i, probing forward when that id is taken. A tax_id shared by N
# tips would otherwise map every trait_df row of those tips onto a single tip,
# producing duplicate `species` rows (with possibly conflicting trait values)
# that break the PSS computation and the CI heatmap.
resolve_duplicate_taxids <- function(df) {
  dup_taxids <- df$tax_id[duplicated(df$tax_id)]
  if (length(dup_taxids) == 0) return(df)

  all_existing <- as.character(unique(df$tax_id))
  resolved <- df

  for (tid in unique(dup_taxids)) {
    species_list <- sort(df$species[df$tax_id == tid])
    kept <- species_list[1]
    duplicates <- species_list[-1]

    debug_log("TAXONOMY CONFLICT: tax_id %s shared by: %s. Assigning synthetic tax_ids to duplicates.",
              tid, paste(species_list, collapse = ", "))

    for (i in seq_along(duplicates)) {
      dup_sp <- duplicates[i]
      synthetic_tid <- as.character(as.integer(tid) + i)
      attempts <- 0
      while (synthetic_tid %in% all_existing && attempts < 1000) {
        synthetic_tid <- as.character(as.integer(synthetic_tid) + 1)
        attempts <- attempts + 1
      }
      if (attempts >= 1000) {
        stop(sprintf("Could not find unused synthetic tax_id for %s (tried 1000 IDs)", dup_sp))
      }
      resolved$tax_id[resolved$species == dup_sp] <- synthetic_tid
      all_existing <- c(all_existing, synthetic_tid)
      debug_log("Synthetic tax_id %s assigned to '%s' (original: %s, kept: %s)",
                synthetic_tid, dup_sp, tid, kept)
    }
  }

  resolved
}

if (nzchar(tax_id_file) && file.exists(tax_id_file)) {
  tax_id_df <- read.csv(tax_id_file, sep = sep_char, stringsAsFactors = FALSE) %>%
    dplyr::mutate(across(everything(), ~ if(is.character(.)) trimws(.) else .))
  if ("tax_id" %in% names(tax_id_df) && "species" %in% names(tax_id_df)) {
    tax_id_df <- tax_id_df %>% dplyr::select(tax_id, species) %>% dplyr::distinct()
    tax_id_df$tax_id <- as.character(tax_id_df$tax_id)
    tax_id_df <- resolve_duplicate_taxids(tax_id_df)
    has.TAX_ID <- TRUE
    debug_log("tax_id_df rows = %d, distinct taxa = %d", nrow(tax_id_df), length(unique(tax_id_df$tax_id)))
  } else {
    message("TAX_ID file provided but missing required columns; proceeding without TAX_ID mapping.")
  }
} else {
  message("No TAX_ID file provided.")
  # A tax_id column in trait_df alone sets has.TAX_ID but does not create
  # tax_id_df, which needs an external file mapping tree-side species names
  # (which may differ from trait_df's naming) to tax_id. Every later
  # `if (has.TAX_ID)` block that reads tax_id_df therefore also checks
  # `exists("tax_id_df")`.
  if ("tax_id" %in% names(trait_df)) {
    has.TAX_ID <- TRUE
  }
}
debug_log("has.TAX_ID = %s", has.TAX_ID)

# With a mapping file, attach its tax_id to the trait table by species name.
if (has.TAX_ID && exists("tax_id_df")) {
  trait_df <- merge(trait_df, tax_id_df, by = "species", all.x = TRUE)
  debug_log("trait_df merged with tax_id_df: rows = %d, missing tax_id = %d",
            nrow(trait_df), sum(is.na(trait_df$tax_id)))
}

# Collapse tax_id.x / tax_id.y (created by the merge when both tables have the column) into one tax_id.
if (has.TAX_ID && !"tax_id" %in% names(trait_df)) {
  tax_id_cols <- intersect(c("tax_id.x", "tax_id.y"), names(trait_df))
  if (length(tax_id_cols) > 0) {
    if (length(tax_id_cols) == 1) {
      trait_df$tax_id <- as.character(trait_df[[tax_id_cols[1]]])
    } else {
      trait_df$tax_id <- dplyr::coalesce(
        as.character(trait_df[[tax_id_cols[1]]]),
        as.character(trait_df[[tax_id_cols[2]]])
      )
    }
    trait_df <- trait_df %>% dplyr::select(-dplyr::all_of(tax_id_cols))
    debug_log("normalized tax_id from merged columns, missing tax_id = %d", sum(is.na(trait_df$tax_id)))
  } else {
    message("TAX_ID requested but no tax_id column found; proceeding without TAX_ID mapping.")
    has.TAX_ID <- FALSE
  }
}

# Leave TAX_ID mode when no tax_id column survives the steps above.
if (has.TAX_ID && !"tax_id" %in% names(trait_df)) {
  message("TAX_ID column missing after setup; proceeding without TAX_ID mapping.")
  has.TAX_ID <- FALSE
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
  if (has.TAX_ID && exists("tax_id_df")) {
    # tax_id_df holds tree-side species names, which may differ from the trait
    # names after taxonomic reclassification (e.g. Nycticebus_pygmaeus ->
    # Xanthonycticebus_pygmaeus), so tips are matched to traits through tax_id.
    tree_tip_to_taxid <- tax_id_df %>%
      dplyr::filter(species %in% tree_species) %>%
      dplyr::rename(tree_name = species)
    debug_log("tree_tip_to_taxid rows = %d, missing tax_id = %d",
              nrow(tree_tip_to_taxid), sum(is.na(tree_tip_to_taxid$tax_id)))

    # tax_ids present in both the tree (via tax_id_df) and the traits (trait_df$tax_id)
    common_tax_ids <- intersect(tree_tip_to_taxid$tax_id, trait_df$tax_id)
    debug_log("common_tax_ids = %d", length(common_tax_ids))

    # Keep the tips whose tax_id is in common_tax_ids.
    tips_to_keep <- tree_tip_to_taxid$tree_name[tree_tip_to_taxid$tax_id %in% common_tax_ids]
    pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
    debug_log("pruned_tree tips (TAX_ID) = %d, nodes = %d", length(pruned_tree$tip.label), pruned_tree$Nnode)
  } else {
    # Without tax_ids, match by species name (spaces read as underscores).
    pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, gsub(" ", "_", trait_df$species)))
    debug_log("pruned_tree tips (species match) = %d, nodes = %d", length(pruned_tree$tip.label), pruned_tree$Nnode)
  }

  # Reduce trait_df to the species of the pruned tree. With tax_ids, the species
  # name is replaced by the tree-side name.
  if (has.TAX_ID && exists("tax_id_df")) {
    trait_df_ori <- trait_df # Table before pruning
    tree_tax_map <- tax_id_df %>%
      dplyr::filter(species %in% pruned_tree$tip.label) %>%
      dplyr::transmute(tax_id, tree_species = species) %>%
      dplyr::distinct(tax_id, .keep_all = TRUE)

    if (nrow(tree_tax_map) > 0) {
      trait_df <- trait_df %>%
        dplyr::left_join(tree_tax_map, by = "tax_id") %>%
        dplyr::mutate(species = dplyr::coalesce(tree_species, species)) %>%
        dplyr::filter(species %in% pruned_tree$tip.label) %>%
        dplyr::select(-tree_species)
      debug_log("trait_df remapped to tree/alignment species using tax_id: rows = %d", nrow(trait_df))
    } else {
      trait_df <- trait_df %>%
        dplyr::filter(tax_id %in% common_tax_ids)
      debug_log("trait_df filtered by common_tax_ids only: rows = %d", nrow(trait_df))
    }
  } else {
    trait_df_ori <- trait_df # Table before pruning
    trait_df <- trait_df %>%
      dplyr::filter(gsub(" ", "_", species) %in% pruned_tree$tip.label)
    debug_log("trait_df after species tree filter rows = %d", nrow(trait_df))
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
