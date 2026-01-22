# ----------------------------------------
# Sample Size Trait Handling
# ----------------------------------------

###################################################################################################
####### If there is a sample size trait (e.g. counts in a prevalence trait (cases/sample)), #######
####### specify it here, as it will allow for a more in-depth analysis.                     #######
###################################################################################################

N_trait <- get_arg(args, 8, "")

# Check if sample size trait is provided and valid
has.N <- FALSE
N.sym <- rlang::sym(N_trait)
if (nzchar(N_trait) && N_trait %in% names(trait_df)) {
  has.N <- TRUE
} else {
  message("No valid sample size trait provided; proceeding without it.")
}