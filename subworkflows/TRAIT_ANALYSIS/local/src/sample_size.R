# sample_size.R — Detect the optional sample-size (n) and case-count (c) trait columns.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Requires `trait_df` (defined by commons.R) and, optionally, `params$n_trait` and
# `params$c_trait`. Defines:
#   n_trait, c_trait  column names (empty string when not given)
#   has.n, has.c      TRUE when the named column exists in `trait_df`; the
#                     reports branch on these flags to add count-based analyses
# =============================================================================

# Fallback logger, used only when commons.R has not defined debug_log().
if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# A prevalence trait (cases / sample) can be analyzed with its counts: n_trait is
# the number of individuals sampled, c_trait the number of observed cases.
n_trait <- if (exists("params")) params$n_trait else "" # Column with the number of individuals sampled
c_trait <- if (exists("params")) params$c_trait else "" # Column with the number of observed cases
debug_log("n_trait = %s", ifelse(nzchar(n_trait), n_trait, "<none>"))
debug_log("c_trait = %s", ifelse(nzchar(c_trait), c_trait, "<none>"))

# A flag is TRUE only when the column is named and present in trait_df.
has.n <- FALSE
if (nzchar(n_trait) && n_trait %in% names(trait_df)) {
  has.n <- TRUE
  debug_log("has.n = TRUE, n missing = %d", sum(is.na(trait_df[[n_trait]])))
} else {
  message("No valid count trait provided; proceeding without it.")
  debug_log("has.n = FALSE")
}

has.c <- FALSE
if (nzchar(c_trait) && c_trait %in% names(trait_df)) {
  has.c <- TRUE
  debug_log("has.c = TRUE, c missing = %d", sum(is.na(trait_df[[c_trait]])))
} else {
  message("No valid sample size trait provided; proceeding without it.")
  debug_log("has.c = FALSE")
}

