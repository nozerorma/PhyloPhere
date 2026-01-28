# ----------------------------------------
# Sample Size Trait Handling
# ----------------------------------------

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

###################################################################################################
####### If there is a sample size trait (e.g. counts in a prevalence trait (cases/sample)), #######
####### specify it here, as it will allow for a more in-depth analysis.                     #######
###################################################################################################

p_trait <- get_arg(args, 8, "")
n_trait <- get_arg(args, 9, "")
debug_log("n_trait = %s", ifelse(nzchar(n_trait), n_trait, "<none>"))
debug_log("p_trait = %s", ifelse(nzchar(p_trait), p_trait, "<none>"))

# Check if sample size trait is provided and valid
has.p <- FALSE
if (nzchar(p_trait) && p_trait %in% names(trait_df)) {
  has.p <- TRUE
  debug_log("has.p = TRUE, p missing = %d", sum(is.na(trait_df[[p_trait]])))
} else {
  message("No valid sample size trait provided; proceeding without it.")
  debug_log("has.p = FALSE")
}

has.n <- FALSE
if (nzchar(n_trait) && n_trait %in% names(trait_df)) {
  has.n <- TRUE
  debug_log("has.n = TRUE, n missing = %d", sum(is.na(trait_df[[n_trait]])))
} else {
  message("No valid count trait provided; proceeding without it.")
  debug_log("has.n = FALSE")
}
