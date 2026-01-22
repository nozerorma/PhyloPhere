# ----------------------------------------
# TAX_ID Trait Handling
# ----------------------------------------

# Check if a TAX_ID column is present, or file mapping taxonomy IDs to species names is provided
has.TAX_ID <- FALSE
# Is there a tax_id file provided?
tax_id_file <- get_arg(args, 9, "")
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

# If TAX_ID is present, merge tax_id_df with trait_df
if (has.TAX_ID && exists("tax_id_df")) {
  trait_df <- merge(trait_df, tax_id_df, by = "tax_id", all.x = TRUE)
}

# Map tree species to taxonomy IDs if needed
if (has.TAX_ID) {
  # Create a mapping of species to tax_id
  tree_ids <- trait_df %>%
    dplyr::filter(gsub(" ", "_", species) %in% tree_species) %>%
    dplyr::select(species, tax_id) %>%
    dplyr::distinct()

# Find common taxonomic IDs between the tree and the traits
  common_tax_ids <- intersect(tree_ids$tax_id, trait_df$tax_id)
  
# Prune the tree to only include species with common taxonomic IDs
  pruned_species <- drop.tip(tree, setdiff(tree$tip.label, tree_ids$species[tree_ids$tax_id %in% common_tax_ids]))
}