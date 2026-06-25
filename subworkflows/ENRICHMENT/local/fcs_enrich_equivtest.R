#!/usr/bin/env Rscript
# =============================================================================
# fcs_enrich_equivtest.R — numerical equivalence test for the vectorized null
# =============================================================================
# Proves fcs_null_enrichstat_vectorized() reproduces RERconverge::fastwilcoxGMTall()$stat
# EXACTLY (the slow path it replaces), across: clean data, ties, NAs, and the
# num.g / bkgenes>2 set-size filtering. If this passes, the permulation p.perm
# produced by the vectorized path is identical (modulo float) to the old loop.
#
# Run:  micromamba run -n phylophere Rscript fcs_enrich_equivtest.R
# =============================================================================

suppressPackageStartupMessages(library(RERconverge))
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "fcs_enrich.R"))

set.seed(42)
TOL <- 1e-8
fails <- 0

# Build a synthetic RERconverge-style gmt object (genesets + geneset.names).
make_gmt <- function(genes, n_sets = 12, min_sz = 3, max_sz = 40) {
  gs <- lapply(seq_len(n_sets), function(i) {
    k <- sample(min_sz:max_sz, 1)
    sample(genes, min(k, length(genes)))
  })
  names(gs) <- paste0("SET_", seq_len(n_sets))
  list(genesets = gs, geneset.names = names(gs))
}

# Run fastwilcoxGMTall (truth) and the vectorized null on the SAME vals, compare stat.
check_case <- function(label, vals, gmt, num_g = 10) {
  gmts <- list(testdb = gmt)
  truth <- tryCatch(RERconverge::fastwilcoxGMTall(vals, gmts, outputGeneVals = FALSE,
                                                  num.g = num_g)$testdb,
                    error = function(e) NULL)
  if (is.null(truth) || nrow(truth) == 0) { cat(sprintf("  [skip] %s: no truth sets\n", label)); return(invisible()) }

  # realenrich must list exactly the sets fastwilcoxGMT reported (it drops NA-pval sets).
  realenrich <- list(testdb = data.frame(pval = truth$pval, stat = truth$stat,
                                          row.names = rownames(truth)))
  corStat <- matrix(vals, ncol = 1, dimnames = list(names(vals), NULL))
  vec <- fcs_null_enrichstat_vectorized(corStat, gmts, realenrich, num_g = num_g)$testdb

  common <- intersect(rownames(truth), rownames(vec))
  d <- abs(truth[common, "stat"] - vec[common, 1])
  maxd <- max(d, na.rm = TRUE)
  # also confirm the NA pattern agrees (sets filtered by num.g/bkgenes)
  na_ok <- all(is.na(vec[common, 1]) == is.na(truth[common, "stat"]) |
               (!is.na(vec[common,1]) & !is.na(truth[common,"stat"])))
  ok <- is.finite(maxd) && maxd < TOL && na_ok
  cat(sprintf("  [%s] %-22s  sets=%d  max|Δstat|=%.2e\n",
              ifelse(ok, "PASS", "FAIL"), label, length(common), maxd))
  if (!ok) fails <<- fails + 1
}

genes <- paste0("G", 1:500)

# 1. Clean continuous scores (no ties, no NA)
vals <- setNames(rnorm(length(genes)), genes)
check_case("clean", vals, make_gmt(genes), num_g = 10)

# 2. Heavy ties (zero-floored magnitude ranking: most genes exactly 0)
vals_tie <- setNames(c(rep(0, 400), abs(rnorm(100))), genes)
check_case("zero-floored-ties", vals_tie, make_gmt(genes), num_g = 5)

# 3. NAs present (genes missing in this null column)
vals_na <- vals; vals_na[sample(length(vals_na), 60)] <- NA
check_case("with-NA", vals_na, make_gmt(genes), num_g = 10)

# 4. Small num.g (more sets pass) and large num.g (filters small sets)
check_case("num_g=3", vals, make_gmt(genes), num_g = 3)
check_case("num_g=30", vals, make_gmt(genes), num_g = 30)

# 5. Multi-column null reproduces each column independently (the real use case).
gmt  <- make_gmt(genes)
gmts <- list(testdb = gmt)
N    <- 5
cols <- replicate(N, setNames(rnorm(length(genes)), genes), simplify = FALSE)
# Reference set list from column 1; truth = fastwilcoxGMTall per column, aligned to it.
ref     <- RERconverge::fastwilcoxGMTall(cols[[1]], gmts, outputGeneVals = FALSE, num.g = 10)$testdb
setrows <- rownames(ref)
truth_stats <- vapply(cols, function(v) {
  t <- RERconverge::fastwilcoxGMTall(v, gmts, outputGeneVals = FALSE, num.g = 10)$testdb
  t[setrows, "stat"]
}, numeric(length(setrows)))
rownames(truth_stats) <- setrows
realenrich    <- list(testdb = data.frame(pval = ref$pval, stat = ref$stat, row.names = setrows))
corStat_multi <- do.call(cbind, cols); rownames(corStat_multi) <- genes
vecN <- fcs_null_enrichstat_vectorized(corStat_multi, gmts, realenrich, num_g = 10)$testdb
maxd_multi <- max(abs(truth_stats[setrows, ] - vecN[setrows, ]), na.rm = TRUE)
ok <- is.finite(maxd_multi) && maxd_multi < TOL
cat(sprintf("  [%s] %-22s  cols=%d  max|Δstat|=%.2e\n",
            ifelse(ok, "PASS", "FAIL"), "multi-column", N, maxd_multi))
if (!ok) fails <- fails + 1

cat(sprintf("\n%s — %d failure(s)\n", ifelse(fails == 0, "ALL PASS", "FAILURES"), fails))
quit(status = if (fails == 0) 0 else 1)
