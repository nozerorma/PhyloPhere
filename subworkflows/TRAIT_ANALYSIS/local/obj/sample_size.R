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

n_trait <- get_arg(args, 8, "") # Trait with number of individuals sampled (population size)
c_trait <- get_arg(args, 9, "") # Trait with number of observed cases (e.g., number of diseased individuals)
debug_log("n_trait = %s", ifelse(nzchar(n_trait), n_trait, "<none>"))
debug_log("c_trait = %s", ifelse(nzchar(c_trait), c_trait, "<none>"))

# Check if sample size trait is provided and valid
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

