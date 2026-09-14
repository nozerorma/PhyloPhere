#!/usr/bin/env Rscript
# Tests for the scoring_compute.R §2g scheme-level aggregation fix — run:
#   Rscript test_scoring_compute_scheme_columns.R
#
# scoring_compute.R is a top-to-bottom pipeline script (no CLI-free entry
# point), so this exercises the exact summarise expressions it uses for
# is_conserved_meta_by_scheme / conserved_pair_by_scheme (copied verbatim from
# the §2g `pos_scores` block) against a synthetic tibble, rather than sourcing
# the whole script.
suppressMessages(library(dplyr))

df <- tibble::tibble(
  Gene = c("G1", "G1"),
  Position = c(10L, 10L),
  side = c("top", "top"),
  caap_group = c("US", "GS3"),
  is_conserved_meta = c(TRUE, FALSE),
  conserved_pair = c("A/B", ""),
  all_mrca_posterior = c(0.9, 0.9),
)

.pos_grp_keys <- c("Gene", "Position", "side")

pos_scores <- df %>%
  group_by(across(all_of(.pos_grp_keys))) %>%
  summarise(
    # Order matters: these must run BEFORE `is_conserved_meta`/`conserved_pair`
    # are overwritten by first() below, or summarise()'s sequential evaluation
    # shadows the per-row vector with the already-collapsed scalar (this
    # exact ordering bug was caught by this test during development).
    is_conserved_meta_by_scheme = paste(
      sort(unique(paste0(as.character(caap_group), ":", as.character(is_conserved_meta)))),
      collapse = ","
    ),
    conserved_pair_by_scheme = paste(
      sort(unique(paste0(as.character(caap_group), ":", conserved_pair)[nzchar(conserved_pair)])),
      collapse = ","
    ),
    is_conserved_meta  = first(is_conserved_meta),
    conserved_pair     = first(conserved_pair),
    all_mrca_posterior = first(all_mrca_posterior),
    .groups = "drop"
  )

ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}

# The display pick (`first()`, US-priority row since df is expected pre-sorted
# by scheme_priority) is unchanged.
check(identical(pos_scores$is_conserved_meta, TRUE),
      "is_conserved_meta display pick == first row (US)")

# The new sibling captures BOTH schemes' values -- this is what the fix adds:
# a display-hidden disagreement (US TRUE vs GS3 FALSE) is no longer silently lost.
check(identical(pos_scores$is_conserved_meta_by_scheme, "GS3:FALSE,US:TRUE"),
      "is_conserved_meta_by_scheme captures both schemes")

check(identical(pos_scores$conserved_pair_by_scheme, "US:A/B"),
      "conserved_pair_by_scheme drops the empty GS3 value, keeps US:A/B")

# all_mrca_posterior (scheme-invariant ASR output) must stay untouched (still
# first()), i.e. the fix must not have accidentally altered this column.
check(identical(pos_scores$all_mrca_posterior, 0.9),
      "all_mrca_posterior unaffected (still first())")

if (!ok) quit(status = 1)
cat("\nAll checks passed.\n")
