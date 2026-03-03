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
# File: rer_objects.R
#

# Set up variable to control command line arguments
args <- commandArgs(TRUE)

# Load libraries
library(dplyr)
library(ape)
library(RERconverge)

# Reading in gene trees with `readTrees`
geneTreesPath <- args[1]
print(geneTreesPath)

# Load trees into gene_trees object
gene_trees <- read.tree(geneTreesPath)

# Prune trees for species in trait file
# Load traitfile
## File must have a multi-column structure with at least a "species" and "trait" column
ori_traits <- read.csv(args[2])
sp_colname <- args[3]

print(sp_colname)

# ── Taxid-based tip-name translation (optional) ──────────────────────────────
# When a taxid mapping file is provided as args[6], rename gene-tree tip labels
# to match the canonical species names used in the traitfile, resolving synonyms
# and name variants via NCBI taxid as a shared key.  This mirrors the same logic
# used for the species tree in TRAIT_ANALYSIS/local/obj/phylo.R.
tax_id_file <- if (length(args) >= 6) args[6] else ""
has_taxid   <- nzchar(tax_id_file) && file.exists(tax_id_file)

if (has_taxid) {
  message("[RER_TREES] Loading taxid file: ", tax_id_file)
  sep_char <- if (endsWith(tax_id_file, ".tsv")) "\t" else ","
  tax_id_df <- read.csv(tax_id_file, sep = sep_char, stringsAsFactors = FALSE)
  tax_id_df[] <- lapply(tax_id_df, function(col) if (is.character(col)) trimws(col) else col)

  # Keep only the two columns we need (tax_id and species)
  if (!all(c("tax_id", "species") %in% names(tax_id_df))) {
    warning("[RER_TREES] taxid file missing 'tax_id' or 'species' column – skipping translation")
    has_taxid <- FALSE
  }
}

if (has_taxid) {
  tax_id_df <- tax_id_df[, c("tax_id", "species"), drop = FALSE]
  tax_id_df <- unique(tax_id_df)

  # Build lookup: any species name (synonym or scientific) -> taxid
  tip_to_taxid <- setNames(tax_id_df$tax_id, gsub(" ", "_", tax_id_df$species))

  # Build canonical name for each taxid: prefer scientific name if that column
  # exists, otherwise take the first entry per taxid.
  # 'canonical' maps taxid -> the species name as used in the traitfile
  trait_species <- unique(gsub(" ", "_", ori_traits[[sp_colname]]))
  # From the taxid lookup, which taxids correspond to traitfile species?
  trait_taxids  <- tip_to_taxid[names(tip_to_taxid) %in% trait_species]
  # Reverse: for each taxid, what is the traitfile canonical name?
  # (use the traitfile name itself, not the taxid file name)
  taxid_to_canonical <- tapply(
    names(trait_taxids), trait_taxids,
    FUN = function(nms) {
      # Prefer the exact traitfile name; if multiple, pick first
      hit <- nms[nms %in% trait_species]
      if (length(hit) > 0) hit[1] else nms[1]
    }
  )
  message(sprintf("[RER_TREES] taxid map: %d unique taxids with traitfile coverage",
                  length(taxid_to_canonical)))
}

# Prune (and optionally rename) the trees

gene_trees_pruned <- list()

for (tree_name in names(gene_trees)) {
    phylotree <- gene_trees[[tree_name]]

    # ── Taxid translation: rename tips before pruning ──────────────────────
    if (has_taxid) {
      old_labels <- phylotree$tip.label
      new_labels <- old_labels
      for (i in seq_along(old_labels)) {
        tip <- old_labels[i]
        if (tip %in% names(tip_to_taxid)) {
          txid <- tip_to_taxid[[tip]]
          canonical <- taxid_to_canonical[[as.character(txid)]]
          if (!is.null(canonical) && !is.na(canonical)) {
            new_labels[i] <- canonical
          }
        }
      }
      phylotree$tip.label <- new_labels
    }
    # ───────────────────────────────────────────────────────────────────────

    gene_trees_pruned[[tree_name]] <- drop.tip(phylotree, setdiff(phylotree$tip.label, ori_traits[[sp_colname]]))
}

print("Prunned")

# Export

pruned_out <- args[4]

# Open a connection to the output file
file_conn <- file(pruned_out, open = "wt")

for (genetree in names(gene_trees_pruned)) {
  # Get the tree
  tree <- gene_trees_pruned[[genetree]]

  # Get the Newick format of the tree as a string
  newick_tree <- write.tree(tree, file = "", append = FALSE)

  # Write the gene name and tree to the output file in a single line without adding a newline
  write(paste0(genetree, "\t", newick_tree), file = file_conn, append = TRUE)
}

# Close the file connection
close(file_conn)

# Pruned trees path
prunedTreesPath <- args[4]
print(prunedTreesPath)

## Generating masterTree
geneTrees <- readTrees(prunedTreesPath) # Set max.read = 100 to toy-up

## Write our masterTree to path
saveRDS(geneTrees, args[5])

### DONE ###
