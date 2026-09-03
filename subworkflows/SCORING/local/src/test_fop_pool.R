#!/usr/bin/env Rscript
# Tests for fop_pool.R — run: Rscript test_fop_pool.R
suppressPackageStartupMessages(library(dplyr))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fop_pool.R"))

ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}
approx <- function(a, b, tol = 1e-6) is.finite(a) && is.finite(b) && abs(a - b) < tol

# ── _p_at_least_2 matches the Python algebra ────────────────────────────────
check(approx(.p_at_least_2(c(1, 1)), 1),                 "_p_at_least_2(1,1) == 1")
check(approx(.p_at_least_2(c(0.5, 0.5)), 0.25),          "_p_at_least_2(.5,.5) == .25")
check(approx(.p_at_least_2(c(0.9, 0.1)), 0.09),          "_p_at_least_2(.9,.1) == .09")
check(approx(.p_at_least_2(c(0.8, 0.8)), 0.64),          "_p_at_least_2(.8,.8) == .64")
check(approx(.p_at_least_2(c(1)), 0),                    "_p_at_least_2 single == 0")
check(approx(.p_at_least_2(c(0.8, 0.8, 0.8)), 0.896),    "_p_at_least_2(.8x3) == .896")

# ── .wmean ─────────────────────────────────────────────────────────────────
check(approx(.wmean(c(1, 0), c(3, 1)), 0.75),            ".wmean weighted")
check(approx(.wmean(c(1, 0), c(NA, NA)), 0.5),           ".wmean falls back to plain mean")
check(approx(.wmean(c(0.4, NA), c(1, 1)), 0.4),          ".wmean drops non-finite x")

# ── single-hypothesis input is a no-op ─────────────────────────────────────
one <- data.frame(
  Gene = "G", Position = 10L, caap_group = "US", hyp_id = NA_character_,
  asr_path_score = 0.42, independence = 0.9, mrca_diversity = 0.5,
  derived_agreement = 1, conservation_gate = 1, core = 0.6,
  mrca_1_path_score = 0.8, mrca_1_node = "n1",
  mrca_2_path_score = 0.7, mrca_2_node = "n2",
  stringsAsFactors = FALSE
)
res1 <- apply_fop_pooling(one, NULL)
check(nrow(res1) == 1, "single-hyp: one row out")
check(approx(res1$asr_path_score, 0.42), "single-hyp: asr_path_score untouched")
check("n_hypotheses" %in% names(res1) && res1$n_hypotheses == 0, "single-hyp: n_hypotheses == 0")

# ── FOP: two hypotheses, two domains ───────────────────────────────────────
# Domain 1: H1 pair (node a, s=0.9), H2 pair (node b, s=0.3)
# Domain 2: H1 pair (node c, s=0.8), H2 pair (node c, s=0.8) [same pair]
fop <- data.frame(
  Gene = "G", Position = 20L, caap_group = "US",
  hyp_id = c("H1", "H2"),
  asr_path_score = c(0.7, 0.2),
  independence = c(1, 1), mrca_diversity = c(1, 0.5),
  derived_agreement = c(1, 1), conservation_gate = c(1, 1),
  core = c(0.72, 0.24),
  mrca_1_path_score = c(0.9, 0.3), mrca_1_node = c("a", "b"),
  mrca_2_path_score = c(0.8, 0.8), mrca_2_node = c("c", "c"),
  stringsAsFactors = FALSE
)
hp <- data.frame(
  hypothesis_id = c("H1", "H1", "H2", "H2"),
  pair = c(1L, 2L, 1L, 2L),
  pss_score = c(10, 10, 2, 2),   # H1 weight 10, H2 weight 2
  stringsAsFactors = FALSE
)
hp_path <- tempfile(fileext = ".tsv")
write.table(hp, hp_path, sep = "\t", quote = FALSE, row.names = FALSE)

resF <- apply_fop_pooling(fop, hp_path)
check(nrow(resF) == 1, "FOP: collapses to one row")
check(resF$n_hypotheses == 2, "FOP: n_hypotheses == 2")
# c_1 = wmean(s=(0.9,0.3), w=(10,2)) = (9 + 0.6)/12 = 0.8
check(approx(resF$pooled_domain_1, 0.8), "FOP: c_1 PSS-weighted == 0.8")
# c_2 = single distinct pair (node c) -> 0.8
check(approx(resF$pooled_domain_2, 0.8), "FOP: c_2 (one distinct pair) == 0.8")
# core = _p_at_least_2(0.8, 0.8) = 0.64
check(approx(resF$core, 0.64), "FOP: pooled core == 0.64")
# diversity pooled = wmean(1, 0.5; w=10,2) = (10 + 1)/12 = 0.91666
check(approx(resF$mrca_diversity, 11/12), "FOP: pooled diversity == 11/12")
# asr = indep(1) * core(0.64) * (0.75 + 0.25*11/12) * da(1) * cg(1)
exp_asr <- 1 * 0.64 * (0.75 + 0.25 * (11/12)) * 1 * 1
check(approx(resF$asr_path_score, exp_asr), sprintf("FOP: pooled asr == %.5f", exp_asr))
# PSS-weighted mean rewards CONSISTENTLY clean domains: H2's weak dom-1 pair
# (kept by the harvest, so real fragility evidence) pulls c_1 below H1's 0.9,
# and the pooled score lands between the two hypotheses rather than at the max.
check(resF$asr_path_score > min(fop$asr_path_score) &&
      resF$asr_path_score < max(fop$asr_path_score),
      "FOP: pooled asr sits between the hypotheses (weighted-mean, not carrier-max)")

# ── FOP without a PSS file -> equal weights ────────────────────────────────
resE <- apply_fop_pooling(fop, NULL)
# c_1 = mean(0.9, 0.3) = 0.6 ; c_2 = 0.8
check(approx(resE$pooled_domain_1, 0.6), "FOP no-PSS: c_1 equal-weight == 0.6")
check(approx(resE$core, .p_at_least_2(c(0.6, 0.8))), "FOP no-PSS: core from equal-weight c")

cat("\n", if (ok) "ALL TESTS PASSED" else "SOME TESTS FAILED", "\n", sep = "")
quit(status = if (ok) 0 else 1)
