# Species tree files and RER gene-tree file of the Tier 2 fixture, from the RERconverge trees object.
#
# usage: Rscript --vanilla rds_to_inputs.R <trees.rds> <species_map.tsv> <outdir>
#   species_map.tsv: columns rer_code (tip label of the object) and species (Genus_species)
# writes in <outdir>:
#   tree_substitution.nwk   master tree of the object, placentals only, rooted, substitution units
#   tree.nwk                the same tree as an ultrametric chronogram (the pipeline assumes a time tree)
#   gene_trees.nwk          one "GENE<TAB>newick" line per gene tree, pruned to the placentals
#   rer_genes.tsv           gene, number of placental species in its tree
suppressPackageStartupMessages(library(ape))
args <- commandArgs(TRUE)
rds <- args[1]; map_file <- args[2]; outdir <- args[3]

map <- read.delim(map_file, stringsAsFactors = FALSE)
nm <- setNames(map$species, map$rer_code)
x <- readRDS(rds)
stopifnot(all(map$rer_code %in% x$masterTree$tip.label))

prune_rename <- function(t) {
  drop <- setdiff(t$tip.label, map$rer_code)
  if (length(drop)) t <- drop.tip(t, drop)
  t$tip.label <- unname(nm[t$tip.label])
  t
}

# Species tree. The master tree is unrooted; the placental root is the split between Atlantogenata
# (Xenarthra + Afrotheria) and the remaining placentals (Meredith et al. 2011 topology).
mt <- prune_rename(x$masterTree)
atlantogenata <- c("Dasypus_novemcinctus", "Orycteropus_afer", "Echinops_telfairi", "Chrysochloris_asiatica",
                   "Elephantulus_edwardii", "Trichechus_manatus", "Loxodonta_africana")
stopifnot(is.monophyletic(mt, atlantogenata))
r <- root(mt, outgroup = atlantogenata, resolve.root = TRUE)
r <- multi2di(r)
r$edge.length[r$edge.length <= 0] <- 1e-8
write.tree(r, file.path(outdir, "tree_substitution.nwk"))
# Same dating as the Tier 1 fixtures: penalised likelihood, lambda = 1, correlated rates, root age 1.
# The input is the tree as read back from the file: the optimizer is sensitive to the internal edge order,
# and the file is what the fixture holds.
r <- read.tree(file.path(outdir, "tree_substitution.nwk"))
u <- suppressWarnings(chronos(r, lambda = 1, model = "correlated", quiet = TRUE,
                              control = chronos.control(eval.max = 10000, iter.max = 10000)))
cat("chronos:", attr(u, "message"), "\n")
stopifnot(!grepl("without", attr(u, "message")))  # a hit of the evaluation limit is not a result
u <- ladderize(structure(unclass(u), class = "phylo"))
stopifnot(is.rooted(u), is.ultrametric(u, tol = 1e-6), Ntip(u) == nrow(map))
write.tree(u, file.path(outdir, "tree.nwk"))

# Gene trees: same pruning and renaming, one line per gene.
genes <- names(x$trees)
n_sp <- integer(length(genes))
con <- file(file.path(outdir, "gene_trees.nwk"), "w")
for (i in seq_along(genes)) {
  t <- prune_rename(x$trees[[i]])
  n_sp[i] <- length(t$tip.label)
  if (n_sp[i] >= 3) writeLines(paste0(genes[i], "\t", write.tree(t)), con)
}
close(con)
writeLines(paste(genes, n_sp, sep = "\t"), file.path(outdir, "rer_genes.tsv"))
cat("tree tips:", Ntip(u), " gene trees:", length(genes), " median placentals per gene:", median(n_sp), "\n")
