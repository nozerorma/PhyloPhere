#!/usr/bin/env Rscript
# test_null_universe.R: the N and the positional universe helpers of scoring_compute.R.
# Run: Rscript tests/null_universe/test_null_universe.R   (scratch files go to ~/.cc/null_universe_test_r, removed afterwards)
suppressPackageStartupMessages({ library(dplyr); library(testthat) })

here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", a[1])))
})
src <- readLines(file.path(here, "..", "..", "subworkflows", "SCORING", "local", "src", "scoring_compute.R"))
i0 <- grep("^# ── Null universe helpers", src)
i1 <- grep("^# ── end of the null universe helpers", src)
stopifnot(length(i0) == 1, length(i1) == 1, i0 < i1)
eval(parse(text = src[i0:i1]))   # the helpers only

work <- file.path(Sys.getenv("HOME"), ".cc", "null_universe_test_r")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE))

test_that(".resolve_n_cycles takes the roster, then caas_perms, then the cycles present", {
  present <- c("b_1", "b_2", "b_3")
  r <- .resolve_n_cycles(present, roster = paste0("b_", 1:5), perms_roster = paste0("b_", 1:4))
  expect_equal(length(r$cycles), 5); expect_equal(r$source, "cycle_roster.txt")
  r <- .resolve_n_cycles(present, roster = NULL, perms_roster = paste0("b_", 1:4))
  expect_equal(length(r$cycles), 4); expect_equal(r$source, "caas_perms.rds")
  msg <- capture.output(r <- .resolve_n_cycles(present, roster = NULL, perms_roster = NULL), type = "message")
  expect_equal(length(r$cycles), 3); expect_match(r$source, "cycles present")
  expect_true(any(grepl("understates N", msg)))          # the fallback is logged
  # a roster that misses a present cycle is skipped, the next source is tried
  msg <- capture.output(r <- .resolve_n_cycles(present, roster = c("b_1", "b_2"), perms_roster = paste0("b_", 1:4)),
                        type = "message")
  expect_equal(r$source, "caas_perms.rds"); expect_true(any(grepl("does not list 1 cycle", msg)))
})

test_that(".cycle_order sorts by numeric suffix", {
  expect_equal(.cycle_order(c("b_10", "b_2", "b_1")), c("b_1", "b_2", "b_10"))
})

test_that(".tested_membership parses only the genes asked for and handles odd lines", {
  f <- file.path(work, "background.output")
  writeLines(c("Gene\tPosition", "G1\t1,2,3,10", "G2\tNULL", "G3\t5, 6", "G4\t"), f)
  m <- .tested_membership(f, c("G1", "G1", "G2", "G3", "G3", "G9"), c(2L, 4L, 1L, 6L, 7L, 1L))
  expect_equal(m$in_tested, c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(m$gene_covered, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))   # G9 has no line
})

test_that("the universe keeps unseen observed positions and drops untested null positions", {
  cyc <- tibble::tribble(~Gene, ~Position, ~cycle,
                         "G1", 1L, "b_1", "G1", 1L, "b_2", "G1", 2L, "b_2",   # tested
                         "G1", 99L, "b_1",                                    # null-positive, not tested
                         "G2", 5L, "b_3")                                     # tested
  null_pos <- distinct(cyc, Gene, Position)
  tested_ok <- c(TRUE, TRUE, FALSE, TRUE)
  obs <- tibble::tibble(Gene = c("G1", "G3"), Position = c(1L, 7L))       # G3:7 is never scored by the null
  out <- .restrict_null_to_universe(cyc, tested_ok, null_pos, obs)
  expect_equal(nrow(distinct(out, Gene, Position)), 3)                    # G1:1, G1:2, G2:5
  expect_false(any(out$Position == 99L))
  # an observed position outside the tested set keeps its null rows
  obs2 <- tibble::tibble(Gene = "G1", Position = 99L)
  out2 <- .restrict_null_to_universe(cyc, tested_ok, null_pos, obs2)
  expect_true(any(out2$Position == 99L))
})

test_that(".null_universe_summary reports N, |U|, unseen and the saturation curve", {
  cyc <- tibble::tribble(~Gene, ~Position, ~cycle,
                         "G1", 1L, "b_1", "G1", 2L, "b_2", "G1", 1L, "b_2", "G2", 3L, "b_10")
  cycles <- paste0("b_", 1:10)
  s <- .null_universe_summary(cyc, cycles, n_obs = 4L, n_unseen = 2L, n_null_raw = 6L, n_cycles_present = 3L,
                              n_source = "cycle_roster.txt", universe_note = "note")
  v <- function(m) s$value[s$metric == m]
  expect_equal(v("N_cycles"), 10)
  expect_equal(v("N_cycles_present_in_null"), 3)
  expect_equal(v("null_positive_positions_in_universe"), 3)
  expect_equal(v("observed_unseen_positions"), 2)
  expect_equal(v("fraction_unseen"), 0.5)
  expect_equal(v("universe_size"), 3 + 2)
  curve <- s[s$metric == "null_positive_positions_cumulative", ]
  expect_equal(curve$cycles, 1:10)                 # 10 checkpoints for N = 10
  expect_equal(curve$value, c(1, 2, 2, 2, 2, 2, 2, 2, 2, 3))   # G1:1 at b_1, G1:2 at b_2, G2:3 at b_10
  # N below 10 gives one checkpoint per cycle, none duplicated
  s2 <- .null_universe_summary(cyc, c("b_1", "b_2", "b_10"), 1L, 0L, 3L, 3L, "x", "y")
  expect_equal(s2$cycles[s2$metric == "null_positive_positions_cumulative"], 1:3)
})
cat("all null-universe helper tests passed\n")
