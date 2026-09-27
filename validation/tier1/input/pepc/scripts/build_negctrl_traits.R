#!/usr/bin/env Rscript
# =============================================================================
# Negative-control traits for the PEPC fixture.
#
# Draws K labelings from the same generator the CAAS permulation null uses
# (subworkflows/CT/local/scripts/permulations.R): phyloq model fit + AIC
# selection on the observed trait, simulation on the (OU-rescaled if selected)
# tree, rank-matching onto the observed values, and Tier-1 Dunn acceptance by
# the lean contrast selector. Each accepted labeling is written as one trait
# column (nc01..ncK). Run through the OBSERVED path of the pipeline, these test
# whether the observed statistic behaves like one more null draw.
#
# Usage:
#   build_negctrl_traits.R <tree.nwk> <my_traits.tsv> <trait_col> <target_pairs>
#                          <pss_top_pct> <K> <seed> <out_traits.tsv> <out_summary.tsv>
# =============================================================================
suppressPackageStartupMessages({ library(ape); library(geiger) })

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)
tree_path    <- args[1]
traits_path  <- args[2]
trait_col    <- args[3]
target_pairs <- as.integer(args[4])
pss_top_pct  <- as.numeric(args[5])
K            <- as.integer(args[6])
seed         <- as.integer(args[7])
out_traits   <- args[8]
out_summary  <- args[9]

script_dir <- {
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  dirname(normalizePath(sub("^--file=", "", hit[1])))
}
ct_scripts <- normalizePath(file.path(script_dir, "../../../../../subworkflows/CT/local/scripts"))
source(file.path(ct_scripts, "pss_core.R"))
source(file.path(ct_scripts, "lean_contrast_selector.R"))

# Same primitive as permulations.R::simulatevec.
simulatevec <- function(namedvec, tree) {
  rm   <- ratematrix(tree, namedvec)
  sims <- sim.char(tree, rm, nsim = 1)
  setNames(as.data.frame(sims)[, 1], rownames(sims))
}

tr <- read.delim(traits_path, stringsAsFactors = FALSE)
vals <- setNames(suppressWarnings(as.numeric(tr[[trait_col]])), tr$species)
vals <- vals[is.finite(vals)]

tree <- read.tree(tree_path)
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, names(vals)))
pruned.tree <- multi2di(pruned.tree, random = FALSE)
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-8
starting.values <- vals[pruned.tree$tip.label]
D <- cophenetic(pruned.tree)

obs_fits       <- fit_models(pruned.tree, starting.values)
selected_model <- select_model(obs_fits, force_model = NULL)
obs_cov        <- covariances_from_fits(pruned.tree, obs_fits)
if (selected_model == "OU") {
  simulation_tree <- rescale(pruned.tree, "OU", as.numeric(obs_fits$OU$opt$alpha))
  simulation_tree$edge.length[simulation_tree$edge.length <= 0] <- 1e-8
} else {
  simulation_tree <- pruned.tree
}
cat(sprintf("model=%s tips=%d fg=%d target_pairs=%d\n", selected_model,
            length(starting.values), sum(starting.values == 1), target_pairs))

set.seed(seed)
rec_value <- unname(sort(starting.values))
obs_fg <- names(starting.values)[starting.values == 1]
out <- list(); draws <- 0L
while (length(out) < K) {
  draws <- draws + 1L
  sim_v   <- simulatevec(starting.values, simulation_tree)
  sim_ord <- order(sim_v)
  pvec    <- setNames(rec_value[order(sim_ord)], names(sim_v))
  e <- tryCatch(
    evaluate_lean_contrast_selection(
      trait_vec = pvec, D = D, target_pairs = target_pairs, tree = pruned.tree,
      cov_bm = obs_cov$BM, cov_ou = obs_cov$OU, selected_model = selected_model,
      top_pct = pss_top_pct, ordinal = TRUE),
    error = function(err) list(tier = 0L))
  if (e$tier == 1L) out[[length(out) + 1L]] <- pvec
}
cat(sprintf("accepted %d Tier-1 labelings in %d draws\n", K, draws))

ids <- sprintf("nc%02d", seq_len(K))
tab <- data.frame(species = tr$species, stringsAsFactors = FALSE)
for (i in seq_len(K)) tab[[ids[i]]] <- unname(out[[i]][tab$species])
tab$family <- tr$family
write.table(tab, out_traits, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Overlap of each control's foreground with the observed trait's foreground.
summ <- do.call(rbind, lapply(seq_len(K), function(i) {
  fg <- names(out[[i]])[out[[i]] == 1]
  x <- out[[i]][names(starting.values)]; y <- starting.values
  data.frame(trait = ids[i], n_fg = length(fg),
             fg_shared_with_obs = length(intersect(fg, obs_fg)),
             jaccard = length(intersect(fg, obs_fg)) / length(union(fg, obs_fg)),
             phi = suppressWarnings(cor(x, y)))
}))
write.table(summ, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)
