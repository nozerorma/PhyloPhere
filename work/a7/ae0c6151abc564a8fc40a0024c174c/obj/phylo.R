# ----------------------------------------
# TAX_ID Trait Handling
# ----------------------------------------

library(ape)

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Check if a TAX_ID column is present, or file mapping taxonomy IDs to species names is provided
has.TAX_ID <- FALSE
# Is there a tax_id file provided?
tax_id_file <- get_arg(args, 10, "")
debug_log("tax_id_file = %s", ifelse(nzchar(tax_id_file), tax_id_file, "<none>"))
# Check if the file is a csv or tsv
if (endsWith(tax_id_file, ".csv")) {
  sep_char <- ","
} else if (endsWith(tax_id_file, ".tsv")) {
  sep_char <- "\t"
} else {
  sep_char <- ","
}

if (nzchar(tax_id_file) && file.exists(tax_id_file)) {
  tax_id_df <- read.csv(tax_id_file, sep = sep_char, stringsAsFactors = FALSE) %>%
    dplyr::mutate(across(everything(), ~ if(is.character(.)) trimws(.) else .))
  if ("tax_id" %in% names(tax_id_df) && "species" %in% names(tax_id_df)) {
    tax_id_df <- tax_id_df %>% dplyr::select(tax_id, species) %>% dplyr::distinct()
    has.TAX_ID <- TRUE
    debug_log("tax_id_df rows = %d, distinct taxa = %d", nrow(tax_id_df), length(unique(tax_id_df$tax_id)))
  } else {
    message("TAX_ID file provided but missing required columns; proceeding without TAX_ID mapping.")
  }
} else {
  message("No TAX_ID file provided.")
  # Check if TAX_ID column is in trait_df
  if ("tax_id" %in% names(trait_df)) {
    has.TAX_ID <- TRUE
  }
}
debug_log("has.TAX_ID = %s", has.TAX_ID)

# If TAX_ID is present, merge tax_id_df with trait_df
if (has.TAX_ID && exists("tax_id_df")) {
  trait_df <- merge(trait_df, tax_id_df, by = "species", all.x = TRUE)
  debug_log("trait_df merged with tax_id_df: rows = %d, missing tax_id = %d",
            nrow(trait_df), sum(is.na(trait_df$tax_id)))
}

# Normalize tax_id column after merge (merge can create tax_id.x / tax_id.y)
if (has.TAX_ID && !"tax_id" %in% names(trait_df)) {
  tax_id_cols <- intersect(c("tax_id.x", "tax_id.y"), names(trait_df))
  if (length(tax_id_cols) > 0) {
    if (length(tax_id_cols) == 1) {
      trait_df$tax_id <- trait_df[[tax_id_cols[1]]]
    } else {
      trait_df$tax_id <- dplyr::coalesce(trait_df[[tax_id_cols[1]]], trait_df[[tax_id_cols[2]]])
    }
    trait_df <- trait_df %>% dplyr::select(-dplyr::all_of(tax_id_cols))
    debug_log("normalized tax_id from merged columns, missing tax_id = %d", sum(is.na(trait_df$tax_id)))
  } else {
    message("TAX_ID requested but no tax_id column found; proceeding without TAX_ID mapping.")
    has.TAX_ID <- FALSE
  }
}

# Guard: if tax_id is still missing, disable TAX_ID mode
if (has.TAX_ID && !"tax_id" %in% names(trait_df)) {
  message("TAX_ID column missing after setup; proceeding without TAX_ID mapping.")
  has.TAX_ID <- FALSE
}

# Map tree species to taxonomy IDs if needed
if (has.TAX_ID) {
  # Create a mapping of species to tax_id
  tree_ids <- trait_df %>%
    dplyr::filter(gsub(" ", "_", species) %in% tree_species) %>%
    dplyr::select(species, tax_id) %>%
    dplyr::distinct()
  debug_log("tree_ids rows = %d, missing tax_id = %d", nrow(tree_ids), sum(is.na(tree_ids$tax_id)))

# Find common taxonomic IDs between the tree and the traits
  common_tax_ids <- intersect(tree_ids$tax_id, trait_df$tax_id)
  debug_log("common_tax_ids = %d", length(common_tax_ids))
  
# Prune the tree to only include species with common taxonomic IDs
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tree_ids$species[tree_ids$tax_id %in% common_tax_ids]))
  debug_log("pruned_tree tips (TAX_ID) = %d, nodes = %d", length(pruned_tree$tip.label), pruned_tree$Nnode)
} else {
  # Prune the tree to only include species present in the trait data
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, gsub(" ", "_", trait_df$species)))
  debug_log("pruned_tree tips (species match) = %d, nodes = %d", length(pruned_tree$tip.label), pruned_tree$Nnode)
}

# Update trait_df to only include species present in the pruned tree
# NOTE: Should I actually do this? Probably in the dataset exploration I want to see all species even if not in the tree.
if (has.TAX_ID) {
  trait_df <- trait_df %>%
    dplyr::filter(tax_id %in% common_tax_ids)
  debug_log("trait_df after TAX_ID tree filter rows = %d", nrow(trait_df))
} else {
  trait_df <- trait_df %>%
    dplyr::filter(gsub(" ", "_", species) %in% pruned_tree$tip.label)
  debug_log("trait_df after species tree filter rows = %d", nrow(trait_df))
}

################### 

# Now we are going to identify the deepest common node for each taxon
# and prepare data for clade labels
# Diagnostic and fixed version
find_taxon_mrca <- function(tree, df) {
  library(ape)
  library(dplyr)
  
  # Function to trace path from tip to root
  trace_to_root <- function(node, tree) {
    path <- node
    current <- node
    edges <- tree$edge
    
    # Maximum iterations to prevent infinite loops
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
  
  # Check if nodes in df are actually in the tree
  tip_labels <- tree$tip.label
  all_nodes <- c(1:length(tip_labels), unique(tree$edge[,1]))
  
  # Get unique taxa
  unique_taxa <- unique(df$taxa)
  unique_taxa <- unique_taxa[!is.na(unique_taxa)]  # Remove NA taxa
  
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
      # Try using ape's built-in MRCA function if nodes are tips
      # Check if these are tip nodes (species names)
      if("species" %in% colnames(df)) {
        taxon_species <- taxon_data$species
        
        # Try to find tips by species names
        tip_indices <- match(taxon_species, tree$tip.label)
        
        if(all(!is.na(tip_indices)) && length(tip_indices) > 1) {
          # Use ape's getMRCA function
          mrca <- getMRCA(tree, tip_indices)
        } else {
          # Fall back to manual path tracing
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
        # No species column, use node tracing
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
