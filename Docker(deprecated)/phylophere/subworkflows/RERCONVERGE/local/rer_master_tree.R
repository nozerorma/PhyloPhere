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

# Prune the trees

gene_trees_pruned <- list()

for (tree_name in names(gene_trees)) {
    phylotree <- gene_trees[[tree_name]]
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
