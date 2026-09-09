#!/usr/bin/env Rscript
# =============================================================================
# T0 shared-oracle golden: fop_pool.R vs the frozen pool_hypotheses output
# (roadmap tier T0, hazard H8 — the R/Python twins have drifted twice).
#
# Reads the SAME file as test_fop_pool_golden.py
#   ../../../CT_DISAMBIGUATION/local/src/convergence/golden/fop_pool_fixture.json
# converts each record-list scenario into the wide disambiguation frame
# apply_fop_pooling() expects, and asserts the pooled axes match the stored
# `expected` (produced by golden/gen_golden.py from the Python twin) to 1e-9.
#
# Each (hyp, domain) cell is given a DISTINCT mrca_<i>_node so pool_group's
# node-dedup is a no-op — that is the regime in which the observed (R) and
# permulation-null (Python) poolers are defined to agree (see fop_pool.py's
# "Residual R/py difference" note).
#
# Run: Rscript test_fop_pool_golden.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "fop_pool.R"))

fix_path <- normalizePath(file.path(
  here, "..", "..", "..",
  "CT_DISAMBIGUATION/local/src/convergence/golden/fop_pool_fixture.json"),
  mustWork = TRUE)
FIX <- fromJSON(fix_path, simplifyVector = FALSE)

TOL <- 1e-9
ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}
approx <- function(a, b) {
  if (is.null(a) && is.null(b)) return(TRUE)
  if (is.null(a) || is.null(b)) return(FALSE)
  is.finite(a) && is.finite(b) && abs(a - b) < TOL
}

# ── record-list scenario -> wide disambiguation frame ────────────────────────
to_frame <- function(sc) {
  recs <- sc$hyp_records
  n <- length(recs)
  doms <- sort(unique(as.integer(unlist(lapply(recs, function(r)
    names(c(r$pair_scores, r$pair_top_scores, r$pair_bottom_scores)))))))
  df <- data.frame(
    Gene = "G", Position = 1L, caap_group = sc$scheme,
    hyp_id = vapply(recs, function(r) r$hyp, character(1)),
    asr_path_score = vapply(recs, function(r) as.numeric(r$asr_path_score %||% NA), numeric(1)),
    independence = vapply(recs, function(r) as.numeric(r$independence %||% NA), numeric(1)),
    mrca_diversity = vapply(recs, function(r) as.numeric(r$mrca_diversity %||% NA), numeric(1)),
    derived_agreement = vapply(recs, function(r) as.numeric(r$derived_agreement %||% NA), numeric(1)),
    conservation_gate = vapply(recs, function(r) as.numeric(r$conservation_gate %||% NA), numeric(1)),
    core = vapply(recs, function(r) as.numeric(r$core %||% NA), numeric(1)),
    stringsAsFactors = FALSE
  )
  getcell <- function(r, fld, d) {
    v <- r[[fld]][[as.character(d)]]
    if (is.null(v)) NA else v
  }
  for (d in doms) {
    di <- as.integer(d)
    df[[sprintf("mrca_%d_path_score", di)]] <-
      vapply(recs, function(r) as.numeric(getcell(r, "pair_scores", d)), numeric(1))
    # distinct node per (hyp, domain) -> dedup is a no-op
    df[[sprintf("mrca_%d_node", di)]] <-
      paste0("n", di, "_", df$hyp_id)
    tp <- vapply(recs, function(r) as.numeric(getcell(r, "pair_top_scores", d)), numeric(1))
    bp <- vapply(recs, function(r) as.numeric(getcell(r, "pair_bottom_scores", d)), numeric(1))
    if (any(is.finite(tp))) df[[sprintf("mrca_%d_top_path_score", di)]] <- tp
    if (any(is.finite(bp))) df[[sprintf("mrca_%d_bot_path_score", di)]] <- bp
    ta <- vapply(recs, function(r) as.character(getcell(r, "pair_derived_top", d)), character(1))
    ba <- vapply(recs, function(r) as.character(getcell(r, "pair_derived_bot", d)), character(1))
    if (any(!is.na(ta) & nzchar(ta)) || any(!is.na(ba) & nzchar(ba))) {
      df[[sprintf("mrca_%d_anc_aa", di)]] <- "A"
      df[[sprintf("mrca_%d_top_aa", di)]] <- ifelse(is.na(ta), "", ta)
      df[[sprintf("mrca_%d_bot_aa", di)]] <- ifelse(is.na(ba), "", ba)
    }
  }
  # conserved pairs -> conserved_<j>_node / _cons (j over the union of pids)
  cps <- lapply(recs, function(r) r$conserved_pair_scores)
  if (any(vapply(cps, function(x) length(x) > 0, logical(1)))) {
    pids <- sort(unique(unlist(lapply(cps, names))))
    for (j in seq_along(pids)) {
      pid <- pids[j]
      df[[sprintf("conserved_%d_node", j)]] <- vapply(recs, function(r) {
        v <- r$conserved_pair_nodes[[pid]]; if (is.null(v)) NA_character_ else as.character(v)
      }, character(1))
      df[[sprintf("conserved_%d_cons", j)]] <- vapply(recs, function(r) {
        v <- r$conserved_pair_scores[[pid]]; if (is.null(v)) NA_real_ else as.numeric(v)
      }, numeric(1))
    }
  }
  df
}

`%||%` <- function(a, b) if (is.null(a)) b else a

pss_file <- function(sc) {
  if (length(sc$pss) == 0) return(NULL)
  hp <- do.call(rbind, lapply(sc$pss, function(row) data.frame(
    hypothesis_id = row[[1]], pair = as.integer(row[[2]]),
    pss_score = as.numeric(row[[3]]), stringsAsFactors = FALSE)))
  p <- tempfile(fileext = ".tsv")
  write.table(hp, p, sep = "\t", quote = FALSE, row.names = FALSE)
  p
}

KEYS <- c("asr_path_score", "core", "independence", "mrca_diversity",
          "derived_agreement", "conservation_gate")

for (sc in FIX) {
  df <- to_frame(sc)
  res <- apply_fop_pooling(df, pss_file(sc))
  check(nrow(res) == 1, sprintf("%s: one row out", sc$name))
  for (k in KEYS) {
    exp <- sc$expected[[k]]
    got <- if (k %in% names(res)) res[[k]][1] else NULL
    check(approx(got, exp),
          sprintf("%s: %s == %s (got %s)", sc$name, k,
                  format(exp %||% NA), format(got %||% NA)))
  }
}

cat("\n", if (ok) "ALL FOP GOLDEN SCENARIOS MATCH" else "SOME CHECKS FAILED", "\n", sep = "")
quit(status = if (ok) 0 else 1)
