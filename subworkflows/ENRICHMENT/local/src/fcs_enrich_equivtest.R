#!/usr/bin/env Rscript
# =============================================================================
# fcs_enrich_equivtest.R - numerical equivalence test for the vectorized null
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
source(file.path(here, "percentile_flags.R"))
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
check_case <- function(label, vals, gmt, num_g = 10, max_g = 0) {
  gmts <- list(testdb = gmt)
  truth <- tryCatch(RERconverge::fastwilcoxGMTall(vals, gmts, outputGeneVals = FALSE,
                                                  num.g = num_g)$testdb,
                    error = function(e) NULL)
  if (is.null(truth) || nrow(truth) == 0) { cat(sprintf("  [skip] %s: no truth sets\n", label)); return(invisible()) }

  # realenrich must list exactly the sets fastwilcoxGMT reported (it drops NA-pval sets).
  realenrich <- list(testdb = data.frame(pval = truth$pval, stat = truth$stat,
                                          row.names = rownames(truth)))
  corStat <- matrix(vals, ncol = 1, dimnames = list(names(vals), NULL))
  vec <- fcs_null_enrichstat_vectorized(corStat, gmts, realenrich, num_g = num_g, max_g = max_g)$testdb

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

# 6. corStat_byrank per-ranking selection via fcs_run_all (the CAAS perms shape).
#    Proves fcs_run_all routes each ranking to ITS OWN null matrix when the RDS
#    carries corStat_byrank, and that the RER single-corStat path is unaffected.
#    Strategy: same observed rankings, two RDS shapes:
#      • byrank : list(global=Mg, top=Mt, bottom=Mb)  (distinct null matrices)
#      • single : list(corStat=Mg)                    (RER-style, shared)
#    Then: global p.perm must MATCH across shapes (both use Mg); top/bottom must
#    DIFFER (byrank uses Mt/Mb, single reuses Mg) - i.e. selection actually happens.
{
  set.seed(7)
  genes6 <- paste0("G", 1:400)
  gmt6   <- make_gmt(genes6, n_sets = 16, min_sz = 4, max_sz = 50)
  gmts6  <- list(testdb = gmt6)
  Np     <- 200
  # Non-negative (zero-floored) magnitude rankings → alternative auto = "greater".
  mkrank <- function(shift) { v <- pmax(0, rnorm(length(genes6)) + shift); setNames(v, genes6) }
  rankings6 <- list(global = mkrank(0), top = mkrank(0.3), bottom = mkrank(-0.3))
  # Three CLEARLY distinct null matrices so a mis-selection is detectable.
  mkmat <- function(scale, shift) {
    m <- matrix(pmax(0, rnorm(length(genes6) * Np) * scale + shift),
                nrow = length(genes6), dimnames = list(genes6, NULL)); m
  }
  Mg <- mkmat(1.0, 0.0); Mt <- mkmat(2.0, 0.5); Mb <- mkmat(0.5, -0.2)

  f_by <- tempfile(fileext = ".rds"); saveRDS(list(corStat_byrank = list(global = Mg, top = Mt, bottom = Mb)), f_by)
  f_sg <- tempfile(fileext = ".rds"); saveRDS(list(corStat = Mg), f_sg)

  e_by <- suppressMessages(fcs_run_all(rankings6, gmts6, num_g = 10, perms_file = f_by))
  e_sg <- suppressMessages(fcs_run_all(rankings6, gmts6, num_g = 10, perms_file = f_sg))

  pp <- function(df, rk) { s <- df[df$ranking == rk, ]; setNames(s$p.perm, s$pathway) }
  cmp <- function(rk) { a <- pp(e_by, rk); b <- pp(e_sg, rk); k <- intersect(names(a), names(b))
                        list(n = length(k), maxd = if (length(k)) max(abs(a[k] - b[k]), na.rm = TRUE) else NA_real_) }

  cg <- cmp("global")  # both shapes use Mg  → identical
  ct <- cmp("top")     # byrank=Mt vs single=Mg → must differ somewhere
  has_perm_by <- any(!is.na(pp(e_by, "global"))) && any(!is.na(pp(e_by, "top"))) && any(!is.na(pp(e_by, "bottom")))

  ok_global <- is.finite(cg$maxd) && cg$maxd < TOL              # selection agrees when matrices match
  ok_top    <- is.finite(ct$maxd) && ct$maxd > TOL              # distinct matrix → distinct p.perm
  ok <- has_perm_by && ok_global && ok_top
  cat(sprintf("  [%s] %-22s  global Δ=%.2e (==0)  top Δ=%.2e (>0)  allperm=%s\n",
              ifelse(ok, "PASS", "FAIL"), "corStat_byrank", cg$maxd, ct$maxd, has_perm_by))
  if (!ok) fails <- fails + 1
}

# 7. Lachenbruch null (Part 1 Fisher/hypergeometric + Part 2 Wilcoxon-among-
#    positive-scorers), vectorized vs. literal per-(set,column) fisher.test/
#    wilcox.test, on a small synthetic corStat_rk + one gmt.
#    Part 1 is EXACT (phyper == fisher.test's one-sided p, closed form) so it
#    gets the file's normal TOL. Part 2 uses a closed-form normal
#    approximation to the Wilcoxon U-statistic UNCONDITIONALLY, while
#    wilcox.test() itself silently switches to an EXACT permutation p-value
#    for small samples with no ties (n1,n2 < 50 by default) - so a gap vs.
#    wilcox.test() is EXPECTED for small n, not a bug. Large n1/n2 (both
#    methods then use the same normal approximation) still gets the tight
#    tolerance, so a real bug in that regime isn't masked.
{
  set.seed(99)
  genes7 <- paste0("G", 1:120)
  gmt7   <- make_gmt(genes7, n_sets = 6, min_sz = 8, max_sz = 40)
  gmts7  <- list(testdb = gmt7)
  Ncol   <- 8
  corStat7 <- matrix(pmax(0, rnorm(length(genes7) * Ncol)), nrow = length(genes7),
                     dimnames = list(genes7, NULL))

  set_names7 <- names(gmt7$genesets)
  realenrich7 <- list(testdb = data.frame(stat = rep(0, length(set_names7)), row.names = set_names7))

  chi1_vec <- fcs_null_lachenbruch_binary_vectorized(corStat7, gmts7, realenrich7,
                                                     num_g = 8, max_g = 40)$testdb
  chi2_vec <- fcs_null_lachenbruch_magnitude_vectorized(corStat7, gmts7, realenrich7,
                                                        max_g = 40)$testdb

  universe7 <- genes7
  max_d1 <- 0; max_d2_small <- 0
  for (s in set_names7) {
    p_genes <- intersect(gmt7$genesets[[s]], universe7)
    n1 <- length(p_genes)
    if (n1 < 8 || n1 > 40) next
    for (j in seq_len(Ncol)) {
      col <- corStat7[, j]
      hits <- names(col[col > 0])
      k1 <- length(intersect(p_genes, hits)); k0 <- n1 - k1
      b1 <- length(hits) - k1; b0 <- (length(universe7) - n1) - b1
      ft <- fisher.test(matrix(c(k1, k0, b1, b0), 2, byrow = TRUE), alternative = "greater")
      chi1_ref <- qchisq(max(ft$p.value, 1e-15), df = 1, lower.tail = FALSE)
      max_d1 <- max(max_d1, abs(chi1_ref - chi1_vec[s, j]), na.rm = TRUE)

      pos <- col[col > 0]; p_pos <- intersect(p_genes, names(pos))
      if (length(p_pos) >= 2 && (length(pos) - length(p_pos)) >= 2) {
        bg_pos <- setdiff(names(pos), p_pos)
        wt <- wilcox.test(pos[p_pos], pos[bg_pos], alternative = "greater")
        chi2_ref <- qchisq(max(wt$p.value, 1e-15), df = 1, lower.tail = FALSE)
        max_d2_small <- max(max_d2_small, abs(chi2_ref - chi2_vec[s, j]), na.rm = TRUE)
      }
    }
  }
  ok1 <- is.finite(max_d1) && max_d1 < TOL
  cat(sprintf("  [%s] lachenbruch-null-part1  max|Δchi1|=%.2e (exact, TOL=%.0e)\n",
              ifelse(ok1, "PASS", "FAIL"), max_d1, TOL))
  cat(sprintf("  [info] lachenbruch-null-part2 (small-n) max|Δchi2|=%.2e (approx vs wilcox.test's own EXACT branch here -- gap expected, not asserted)\n",
              max_d2_small))
  if (!ok1) fails <- fails + 1

  # Dedicated large-n scenario: wilcox.test() only switches to ITS OWN normal
  # approximation once both samples are >= 50 (no ties) -- genes7's small
  # universe can't reach that on both sides at once, so build one big enough
  # to. Once both sides use the SAME normal-approximation formula (with
  # wilcox.test's continuity correction, which the vectorized Part 2 now also
  # applies), agreement should be tight.
  set.seed(100)
  # All-positive draws (no zero-floor shrinkage) so the realized positive-
  # scorer counts equal n1L/n2L exactly, comfortably clearing the >= 50
  # threshold on both sides.
  n1L <- 70; n2L <- 100
  xL <- abs(rnorm(n1L, mean = 0.3)) + 0.05; yL <- abs(rnorm(n2L, mean = 0)) + 0.05
  genesL <- c(paste0("S", seq_len(n1L)), paste0("B", seq_len(n2L)))
  valsL  <- setNames(c(xL, yL), genesL)
  # A second, throwaway geneset covering the B genes is required so that
  # fcs_null_enrichstat_vectorized's genes_db (the union of ALL genesets in
  # the GMT, matching fastwilcoxGMT's own background convention) includes a
  # background distinct from SETL -- with only one geneset, genes_db would
  # equal SETL itself and n2 would be 0 (NA), which is a test-construction
  # artifact, not a real degenerate case.
  gmtL <- list(genesets = list(SETL = paste0("S", seq_len(n1L)),
                               OTHER = paste0("B", seq_len(n2L))),
              geneset.names = c("SETL", "OTHER"))
  corStatL <- matrix(valsL, ncol = 1, dimnames = list(genesL, NULL))
  realenrichL <- list(testdb = data.frame(stat = 0, row.names = "SETL"))
  chi2L <- fcs_null_lachenbruch_magnitude_vectorized(corStatL, list(testdb = gmtL), realenrichL,
                                                     max_g = n1L)$testdb
  posL <- valsL[valsL > 0]
  pL   <- intersect(paste0("S", seq_len(n1L)), names(posL))
  bgL  <- setdiff(names(posL), pL)
  wtL  <- wilcox.test(posL[pL], posL[bgL], alternative = "greater")
  chi2L_ref <- qchisq(max(wtL$p.value, 1e-15), df = 1, lower.tail = FALSE)
  max_d2_large <- abs(chi2L_ref - chi2L["SETL", 1])
  ok2 <- is.finite(max_d2_large) && max_d2_large < TOL
  cat(sprintf("  [%s] lachenbruch-null-part2 (large-n, n1=%d n2=%d) max|Δchi2|=%.2e (TOL=%.0e)\n",
              ifelse(ok2, "PASS", "FAIL"), length(pL), length(bgL), max_d2_large, TOL))
  if (!ok2) fails <- fails + 1
}

# 8. Path-sum with null_mat= supplied: (a) reproduces M %*% null_mat sums/NES
#    exactly (the downstream math is untouched by this change - this tests
#    the alignment/plumbing in fcs_run_permulation's new branch, not new
#    statistics), and (b) a gene present in `vals` but absent from
#    null_mat's rownames gets a neutral (0) null row and does not perturb
#    results for sets that don't contain it.
{
  set.seed(123)
  genes8 <- paste0("G", 1:200)
  gmt8   <- make_gmt(genes8, n_sets = 8, min_sz = 5, max_sz = 40)
  gmts8  <- list(testdb = gmt8)
  vals8  <- setNames(pmax(0, rnorm(length(genes8))), genes8)
  Ncol8  <- 50
  null8  <- matrix(pmax(0, rnorm(length(genes8) * Ncol8)), nrow = length(genes8),
                   dimnames = list(genes8, NULL))

  res8 <- fcs_run_permulation(vals8, gmts8, num_g = 5, max_g = 40, null_mat = null8)

  # Hand-computed reference using the exact same M %*% null_mat formula.
  gs8 <- gmt8$genesets; names(gs8) <- gmt8$geneset.names
  valid8 <- sapply(gs8, function(g) length(intersect(g, genes8)) >= 5 && length(intersect(g, genes8)) <= 40)
  gs8v <- gs8[valid8]
  ref_rows <- lapply(names(gs8v), function(pname) {
    g_in <- intersect(gs8v[[pname]], genes8)
    obs_sum <- sum(vals8[g_in])
    null_sums <- colSums(null8[g_in, , drop = FALSE])
    p_ref <- (sum(null_sums >= obs_sum) + 1) / (Ncol8 + 1)
    c(pathway = pname, perm_pval = p_ref)
  })
  ref8 <- data.frame(pathway = sapply(ref_rows, `[[`, "pathway"),
                     perm_pval_ref = as.numeric(sapply(ref_rows, `[[`, "perm_pval")))
  cmp8 <- merge(res8, ref8, by = "pathway")
  maxd8 <- max(abs(cmp8$perm_pval - cmp8$perm_pval_ref), na.rm = TRUE)
  n_perms_ok <- all(res8$perm_pval >= 1 / (Ncol8 + 1) - 1e-12)
  ok8 <- is.finite(maxd8) && maxd8 < TOL && nrow(cmp8) == nrow(res8) && n_perms_ok
  cat(sprintf("  [%s] pathsum-null-mat-plumbing  sets=%d  max|Δp|=%.2e  (N=%d from null_mat, not n_perms_sum)\n",
              ifelse(ok8, "PASS", "FAIL"), nrow(cmp8), maxd8, Ncol8))
  if (!ok8) fails <- fails + 1

  # Neutral-row check: a gene in vals8 but absent from null_mat must not
  # change results for sets that don't contain it.
  extra_gene <- "GX_not_in_null"
  vals8b <- c(vals8, setNames(0.7, extra_gene))
  res8b <- fcs_run_permulation(vals8b, gmts8, num_g = 5, max_g = 40, null_mat = null8)
  # No geneset in make_gmt() can contain extra_gene (drawn only from genes8),
  # so every set's result must be identical with/without the extra gene.
  cmp8b <- merge(res8[, c("pathway", "perm_pval", "perm_nes")],
                res8b[, c("pathway", "perm_pval", "perm_nes")],
                by = "pathway", suffixes = c("_a", "_b"))
  maxd8b <- max(abs(cmp8b$perm_pval_a - cmp8b$perm_pval_b),
               abs(cmp8b$perm_nes_a - cmp8b$perm_nes_b), na.rm = TRUE)
  ok8b <- is.finite(maxd8b) && maxd8b < TOL
  cat(sprintf("  [%s] pathsum-null-mat-neutral-row  max|Δ|=%.2e\n", ifelse(ok8b, "PASS", "FAIL"), maxd8b))
  if (!ok8b) fails <- fails + 1
}

cat(sprintf("\n%s - %d failure(s)\n", ifelse(fails == 0, "ALL PASS", "FAILURES"), fails))
quit(status = if (fails == 0) 0 else 1)
