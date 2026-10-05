#!/usr/bin/env Rscript
# The vectorized permulation null of the FCS against RERconverge::fastwilcoxGMTall, the function the observed side runs.
#
#   Rscript validation/unification/check_fcs_null_equivalence.R [repo root]
#
# For every column of a null matrix with the structure of the CAAS null (mostly structural zeros, rounded scores so that
# ties are heavy, optionally NA cells) the statistic of every set is computed twice: column by column through
# fcs_run_ranking's own call (the max_g filter, then fastwilcoxGMTall), and for all columns at once by
# fcs_null_enrichstat_vectorized. Cells must agree in their NA pattern and to 1e-12. Exit status 1 if any does not.
# Needs RERconverge (installed in the pipeline's environment, not on every development machine).
suppressPackageStartupMessages({ library(RERconverge); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
src <- file.path(root, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "percentile_flags.R"))     # fcs_enrich.R skips its cwd-relative lookup when this is already loaded
source(file.path(src, "fcs_enrich.R"))

make_case <- function(seed, with_na) {
  set.seed(seed)
  G <- 400; N <- 40; genes <- sprintf("g%03d", seq_len(G))
  cs <- matrix(0, G, N, dimnames = list(genes, NULL))
  nz <- matrix(runif(G * N) < 0.15, G, N); cs[nz] <- round(runif(sum(nz)), 1)
  if (with_na) cs[matrix(runif(G * N) < 0.03, G, N)] <- NA
  mk <- function(prefix, k, sizes) {
    gs <- lapply(seq_len(k), function(i) sample(genes[1:360], sample(sizes, 1)))
    names(gs) <- sprintf("%s%02d", prefix, seq_len(k)); list(genesets = gs, geneset.names = names(gs))
  }
  list(cs = cs, gmts = list(dbA = mk("A", 25, c(3, 8, 10, 15, 40, 90, 150)), dbB = mk("B", 12, c(10, 30, 60, 200))))
}

observed_side <- function(vals, gmts, num_g, max_g) {
  # fcs_run_ranking: sets over max_g are dropped before the call; fastwilcoxGMT drops NA values itself
  keep_gmts <- lapply(gmts, function(gmt) {
    gs <- gmt$genesets; names(gs) <- gmt$geneset.names
    keep <- vapply(gs, function(set) length(intersect(set, names(vals))) <= max_g, logical(1))
    list(genesets = gs[keep], geneset.names = gmt$geneset.names[keep])
  })
  RERconverge::fastwilcoxGMTall(vals, keep_gmts, outputGeneVals = FALSE, num.g = num_g)
}

fail <- FALSE
for (case in list(list(seed = 11, na = FALSE), list(seed = 12, na = TRUE))) {
  d <- make_case(case$seed, case$na); num_g <- 10; max_g <- 100
  real <- lapply(d$gmts, function(g) data.frame(pval = rep(NA_real_, length(g$geneset.names)), stat = NA_real_, row.names = g$geneset.names))
  vec <- fcs_null_enrichstat_vectorized(d$cs, d$gmts, real, num_g, max_g)
  cells <- 0; both_defined <- 0; na_mismatch <- 0; worst <- 0
  for (j in seq_len(ncol(d$cs))) {
    ref <- observed_side(d$cs[, j], d$gmts, num_g, max_g)
    for (db in names(d$gmts)) {
      r <- ref[[db]]
      for (s in rownames(real[[db]])) {
        cells <- cells + 1
        rv <- if (!is.null(r) && s %in% rownames(r)) r[s, "stat"] else NA_real_
        vv <- vec[[db]][s, j]
        if (is.na(rv) != is.na(vv)) na_mismatch <- na_mismatch + 1
        else if (!is.na(rv)) { both_defined <- both_defined + 1; worst <- max(worst, abs(rv - vv)) }
      }
    }
  }
  ok <- na_mismatch == 0 && both_defined > 200 && worst < 1e-12
  cat(sprintf("seed %d, NA cells %s: %d cells, %d defined in both, %d NA-pattern mismatches, max |difference| %s -> %s\n",
              case$seed, case$na, cells, both_defined, na_mismatch, format(worst, digits = 3), if (ok) "PASS" else "FAIL"))
  fail <- fail || !ok
}
quit(status = as.integer(fail))
