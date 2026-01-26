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

N_trait <- get_arg(args, 8, "")
debug_log("N_trait = %s", ifelse(nzchar(N_trait), N_trait, "<none>"))

# Check if sample size trait is provided and valid
has.N <- FALSE
if (nzchar(N_trait) && N_trait %in% names(trait_df)) {
  has.N <- TRUE
  debug_log("has.N = TRUE, N missing = %d", sum(is.na(trait_df[[N_trait]])))
} else {
  message("No valid sample size trait provided; proceeding without it.")
  debug_log("has.N = FALSE")
}
