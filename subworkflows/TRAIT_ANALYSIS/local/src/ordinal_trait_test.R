#!/usr/bin/env Rscript
# =============================================================================
# ordinal_trait_test.R — regression test for binary / ordinal fg-bg traits in
# contrast selection. Run: Rscript subworkflows/TRAIT_ANALYSIS/local/src/ordinal_trait_test.R
#
# Covers the fix for: a binary 0/1 (or few-level integer) trait was quantile-
# discretised, `quantile([mostly 0s], 0.9) == 0` labelled every species "top",
# leaving 0 candidate pairs and crashing pair_sel.f with "subscript out of bounds".
# =============================================================================

suppressWarnings(suppressMessages({library(ape); library(dplyr); library(tidyr); library(tibble)}))

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "subworkflows/TRAIT_ANALYSIS/local/src"

has.n <- FALSE
contrast_max_iter <- 3L
debug_log <- function(...) invisible(NULL)
suppressWarnings(source(file.path(here, "stats.R")))
suppressWarnings(source(file.path(here, "selection_algorithm.R")))
suppressWarnings(source(file.path(here, "..", "..", "..", "CT", "local", "scripts", "lean_contrast_selector.R")))

fails <- 0L
ok <- function(cond, msg) {
  cat(if (isTRUE(cond)) "  ok   " else "  FAIL ", msg, "\n", sep = "")
  if (!isTRUE(cond)) fails <<- fails + 1L
}

# ── is_ordinal_trait / compute_trait_thresholds ─────────────────────────────
cat("is_ordinal_trait + compute_trait_thresholds\n")
bin_minority <- c(rep(0, 230), rep(1, 7))
ok(is_ordinal_trait(bin_minority), "binary minority auto-detected as ordinal")
ok(!is_ordinal_trait(rbeta(40, 2, 20)), "continuous (beta) not flagged ordinal")
ok(is_ordinal_trait(c(rep(0, 20), rep(1, 10), rep(2, 5))), "3-level integer flagged ordinal")
ok(!is_ordinal_trait(bin_minority, "continuous"), "trait_type='continuous' forces off")
ok(is_ordinal_trait(rbeta(40, 2, 20), "ordinal"), "trait_type='ordinal' forces on")

th <- compute_trait_thresholds(bin_minority, discrete_method = "parameterized",
                               top_quantile = 0.98, bottom_quantile = 0.45)
ok(identical(th$discrete_method, "ordinal"), "binary -> method 'ordinal'")
ok(th$lower_thresh == 0 && th$upper_thresh == 1, "binary thresholds = {0, 1}")
ok(sum(bin_minority >= th$upper_thresh) == 7, "7 foreground (top) species")
ok(sum(bin_minority <= th$lower_thresh) == 230, "230 background (bottom) species")

th3 <- compute_trait_thresholds(c(rep(0, 20), rep(1, 10), rep(2, 5)))
ok(th3$lower_thresh == 0 && th3$upper_thresh == 2, "3-level thresholds = {0, 2} (mid excluded)")

# ── lean_thresholds (permulation path parity) ──────────────────────────────
cat("lean_thresholds (permulation selector)\n")
lt <- lean_thresholds(setNames(bin_minority, paste0("sp", 1:237)),
                      "parameterized", 0.45, 0.98, "auto")
ok(lt$lower == 0 && lt$upper == 1, "lean_thresholds(binary) = {0, 1}")

# ── evaluate_lean_contrast_selection end to end ────────────────────────────
cat("evaluate_lean_contrast_selection\n")
set.seed(42)
tr <- rtree(40); tr$tip.label <- paste0("sp", 1:40)
D <- cophenetic(tr)
fg_idx <- c(3, 11, 19, 27, 35, 40)
v <- setNames(rep(0, 40), tr$tip.label); v[fg_idx] <- 1
e <- evaluate_lean_contrast_selection(v, D, target_pairs = 3, trait_type = "auto")
ok(e$tier %in% c(1L, 2L), "binary trait -> a pool-eligible design (tier 1/2)")
ok(all(v[e$fg] == 1) && all(v[e$bg] == 0), "fg members are all 1, bg members all 0")

v3 <- setNames(rep(1, 40), tr$tip.label); v3[fg_idx] <- 2; v3[c(1, 5, 9, 13, 17, 21)] <- 0
e3 <- evaluate_lean_contrast_selection(v3, D, target_pairs = 3, trait_type = "auto")
ok(all(v3[e3$fg] == 2) && all(v3[e3$bg] == 0), "3-level: fg = top code, bg = bottom code")

# ── pair_sel.f: binary happy path + empty-candidate guard ──────────────────
cat("pair_sel.f\n")
set.seed(1)
tr2 <- rtree(30); tr2$tip.label <- paste0("sp", 1:30)
D2 <- cophenetic(tr2)
val <- rep(0, 30); val[c(4, 10, 16, 22, 28)] <- 1
traits_df <- data.frame(species = tr2$tip.label, trait = "t", value = val, stringsAsFactors = FALSE)
pw <- expand.grid(species1 = tr2$tip.label, species2 = tr2$tip.label, stringsAsFactors = FALSE)
pw <- pw[pw$species1 != pw$species2, ]
pw$trait_diff <- val[match(pw$species1, tr2$tip.label)] - val[match(pw$species2, tr2$tip.label)]
cat1 <- ifelse(val == 1, "top", "bottom"); names(cat1) <- tr2$tip.label
keep <- cat1[pw$species1] != cat1[pw$species2]
ov <- pw[keep, c("species1", "species2", "trait_diff")]
res <- pair_sel.f(D2, ov, traits_df, "t")
ok(nrow(res$selected_pairs) >= 3, "binary trait -> >=3 contrast pairs selected")
res0 <- suppressWarnings(pair_sel.f(D2, ov[0, ], traits_df, "t"))
ok(nrow(res0$selected_pairs) == 0, "empty overlap_df -> 0 pairs, no crash")

cat(sprintf("\n%s (%d failure%s)\n", if (fails == 0) "PASS" else "FAIL",
            fails, if (fails == 1) "" else "s"))
quit(status = if (fails == 0) 0L else 1L)
