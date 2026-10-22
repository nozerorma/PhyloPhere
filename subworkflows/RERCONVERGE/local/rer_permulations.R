# RERconverge phenotype generators with the same correlation settings for observed
# and null paths. Simulation follows PermulationFuncs.R, commit 2bd328f; only the
# correlation and foreground-path settings are supplied explicitly by the caller.

rer_winsor_limit <- function(x) {
  if (is.null(x) || is.na(x) || x <= 0) NULL else as.numeric(x)
}

rer_score_paths <- function(RERmat, paths, mode, min.sp, min.pos = 2,
                            winR = NULL, winT = NULL) {
  if (mode == "continuous") {
    RERconverge::correlateWithContinuousPhenotype(
      RERmat, paths, min.sp = min.sp,
      winsorizeRER = rer_winsor_limit(winR), winsorizetrait = rer_winsor_limit(winT))
  } else if (mode == "binary") {
    RERconverge::correlateWithBinaryPhenotype(
      RERmat, paths, min.sp = min.sp, min.pos = min.pos, weighted = "auto",
      winsorizeRER = rer_winsor_limit(winR), winsorizetrait = NULL)
  } else stop("Unknown RER mode: ", mode)
}

# The output has the corP/corRho/corStat interface of getPermsContinuous/Binary.
# A missing association remains NA within that column; measured zeros remain zero.
rer_null_columns <- function(numperms, genes, score_column) {
  stopifnot(length(numperms) == 1, is.finite(numperms), numperms >= 1,
            numperms == as.integer(numperms))
  p <- rho <- stat <- matrix(NA_real_, length(genes), numperms,
                             dimnames = list(genes, NULL))
  for (j in seq_len(numperms)) {
    result <- score_column(j)
    stopifnot(setequal(rownames(result), genes))
    result <- result[genes, , drop = FALSE]
    p[, j] <- result$P
    rho[, j] <- result$Rho
    stat[, j] <- sign(result$Rho) * -log10(result$P)
  }
  list(corP = as.data.frame(p), corRho = as.data.frame(rho), corStat = as.data.frame(stat))
}

# simpermvec preserves the observed trait distribution and uses Brownian-motion
# simulation to determine its phylogenetic assignment, as getNullCor(type="simperm").
rer_perms_continuous <- function(numperms, traitvec, RERmat, trees, mastertree,
                                 min.sp, winR = 3, winT = 3) {
  rer_null_columns(numperms, rownames(RERmat), function(j) {
    vec <- RERconverge::simpermvec(traitvec, mastertree)
    paths <- RERconverge::char2Paths(stats::na.omit(vec), trees)
    rer_score_paths(RERmat, paths, "continuous", min.sp, winR = winR, winT = winT)
  })
}

# The CC histories are drawn by the native categoricalPermulations implementation.
# The observed and every null history use the caller's clade and correlation options.
rer_perms_binary <- function(numperms, fg_vec, RERmat, trees, min.sp,
                             min.pos = 2, winR = NULL, clade = "all") {
  stopifnot(clade %in% c("all", "terminal", "ancestral"))
  phenotype <- setNames(rep(0, length(trees$masterTree$tip.label)), trees$masterTree$tip.label)
  phenotype[names(phenotype) %in% fg_vec] <- 1
  simulated <- RERconverge::categoricalPermulations(trees, phenotype, rm = "ER", rp = "auto", ntrees = numperms)
  perm_trees <- lapply(simulated$trees, function(x) {
    tr <- trees$masterTree
    tr$edge.length <- c(x$tips, x$nodes)[tr$edge[, 2]] - 1
    names(tr$edge.length) <- NULL
    tr
  })
  # This helper is internal in the pinned package and is also used by getPermsBinary.
  foregrounds <- getFromNamespace("getForegroundsFromBinaryTree", "RERconverge")
  perm_fg <- lapply(perm_trees, foregrounds)
  stopifnot(length(perm_fg) == numperms)
  rer_null_columns(numperms, rownames(RERmat), function(j) {
    paths <- RERconverge::foreground2Paths(perm_fg[[j]], trees, clade = clade)
    rer_score_paths(RERmat, paths, "binary", min.sp, min.pos, winR = winR)
  })
}
