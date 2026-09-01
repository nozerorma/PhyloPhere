################################################################################
# Analytical phylogenetic shift score
#
# VENDORED VERBATIM from Fabio Barteri's phyloq:
#   github.com/linudz/phyloq @ 84eb428
#   fabio/primate.traits/scripts/pss.core.R
#   (+ fit_models / covariances_from_fits from parametric_bootstrap_chunk.R)
#
# Local additions are minimal and marked "PHYLOPHERE:":
#   - `force_model` arg on phylogenetic_shift_score() / select_model() to honour
#     an explicit --perm_strategy (ou|bm) instead of the AIC rule.
# The BM/OU fit, validation (validate_fit), and model selection
# (select_model: OU iff AIC_OU + 2 < AIC_BM, otherwise BM — plain AIC) MUST stay
# identical to upstream so PhyloPhere's PSS ranking matches phyloq's.
################################################################################

# PHYLOPHERE: fit both models + validate (verbatim from parametric_bootstrap_chunk.R)
fit_models <- function(tree, values) {
  fits <- list(
    BM = suppressWarnings(geiger::fitContinuous(tree, values, model = "BM", ncores = 1)),
    OU = suppressWarnings(geiger::fitContinuous(tree, values, model = "OU", ncores = 1))
  )
  validate_fit(fits$BM, "BM")
  validate_fit(fits$OU, "OU")
  fits
}

# PHYLOPHERE: BM + OU covariance list (verbatim from parametric_bootstrap_chunk.R)
covariances_from_fits <- function(tree, fits) {
  list(
    BM = bm_covariance(tree, fits$BM$opt$sigsq),
    OU = ou_covariance(tree, fits$OU$opt$alpha, fits$OU$opt$sigsq)
  )
}

phylogenetic_shift_score <- function(data, tree, trait,
                                      species_col = "SpeciesBROAD",
                                      outdir = NULL,
                                      write_output = !is.null(outdir),
                                      verbose = TRUE,
                                      force_model = NULL) {
  require_package("ape")
  require_package("geiger")
  say <- function(...) if (isTRUE(verbose)) message(...)
  traits <- read_trait_data(data, species_col)
  phy <- read_phylogeny(tree)
  matched <- match_trait_to_tree(traits, phy, trait, species_col)
  traits <- matched$data
  phy <- matched$tree
  values <- stats::setNames(traits[[trait]], traits[[species_col]])
  values <- values[phy$tip.label]

  say("Trait: ", trait, "; matched species: ", length(values))
  say("Fitting BM and OU.")
  fits <- list(
    BM = geiger::fitContinuous(phy, values, model = "BM", ncores = 1),
    OU = geiger::fitContinuous(phy, values, model = "OU", ncores = 1)
  )
  validate_fit(fits$BM, "BM")
  validate_fit(fits$OU, "OU")

  covariance <- list(
    BM = bm_covariance(phy, fits$BM$opt$sigsq),
    OU = ou_covariance(phy, fits$OU$opt$alpha, fits$OU$opt$sigsq)
  )
  aic_bm <- fit_aic(fits$BM)
  aic_ou <- fit_aic(fits$OU)
  delta_aic <- aic_bm - aic_ou
  selected_model <- select_model(fits, force_model = force_model)
  scores <- calculate_pairwise_scores(
    values, phy, covariance$BM, covariance$OU, selected_model
  )
  scores$AIC_BM <- aic_bm
  scores$AIC_OU <- aic_ou
  scores$delta_AIC <- delta_aic
  scores$SelectedModel <- selected_model

  fit_table <- model_fit_table(fits, trait, selected_model)
  result <- list(
    trait = trait, selected_model = selected_model, tree = phy, fits = fits,
    model_fit = fit_table, scores = scores
  )
  class(result) <- "phylogenetic_shift_score_result"

  if (isTRUE(write_output)) {
    if (is.null(outdir)) stop("outdir is required when write_output is TRUE.", call. = FALSE)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    output_file <- file.path(outdir, paste0(trait, ".score_results.tsv"))
    utils::write.table(scores, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    say("AIC BM: ", aic_bm, "; AIC OU: ", aic_ou, "; delta AIC: ", delta_aic)
    say("Selected model: ", selected_model)
    say("Output: ", output_file)
  }
  result
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required.", call. = FALSE)
  }
}

normalize_species_names <- function(x) {
  x <- gsub("'", "", as.character(x), fixed = TRUE)
  gsub(" ", "_", x, fixed = TRUE)
}

read_trait_data <- function(data, species_col) {
  if (is.character(data) && length(data) == 1L) {
    data <- utils::read.table(data, header = TRUE, sep = "\t",
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(data)) stop("data must be a data frame or TSV path.", call. = FALSE)
  if (!species_col %in% names(data)) stop("Species column '", species_col, "' not found.", call. = FALSE)
  data[[species_col]] <- normalize_species_names(data[[species_col]])
  data
}

read_phylogeny <- function(tree) {
  if (inherits(tree, "phylo")) phy <- tree
  else if (is.character(tree) && length(tree) == 1L) phy <- ape::read.tree(tree)
  else stop("tree must be a phylo object or Newick path.", call. = FALSE)
  phy$tip.label <- normalize_species_names(phy$tip.label)
  phy
}

match_trait_to_tree <- function(data, tree, trait, species_col) {
  if (!trait %in% names(data)) stop("Trait '", trait, "' not found.", call. = FALSE)
  trait_data <- data[, c(species_col, trait), drop = FALSE]
  trait_data[[trait]] <- suppressWarnings(as.numeric(trait_data[[trait]]))
  trait_data <- trait_data[!is.na(trait_data[[species_col]]) & is.finite(trait_data[[trait]]), , drop = FALSE]
  if (anyDuplicated(trait_data[[species_col]])) stop("Duplicated species names in trait data.", call. = FALSE)
  common <- intersect(tree$tip.label, trait_data[[species_col]])
  if (length(common) < 3L) stop("Fewer than three matched species.", call. = FALSE)
  to_drop <- setdiff(tree$tip.label, common)
  if (length(to_drop)) tree <- ape::drop.tip(tree, to_drop)
  trait_data <- trait_data[match(tree$tip.label, trait_data[[species_col]]), , drop = FALSE]
  list(data = trait_data, tree = tree)
}

validate_fit <- function(fit, model) {
  if (is.null(fit$opt$sigsq) || !is.finite(fit$opt$sigsq) || fit$opt$sigsq <= 0) {
    stop(model, " produced an invalid variance rate.", call. = FALSE)
  }
  if (model == "OU" && (is.null(fit$opt$alpha) || !is.finite(fit$opt$alpha) || fit$opt$alpha <= 0)) {
    stop("OU produced an invalid alpha.", call. = FALSE)
  }
}

bm_covariance <- function(tree, sigma2) {
  covariance <- sigma2 * ape::vcv.phylo(tree)
  covariance[tree$tip.label, tree$tip.label, drop = FALSE]
}

ou_covariance <- function(tree, alpha, sigma2) {
  shared <- ape::vcv.phylo(tree)
  shared <- shared[tree$tip.label, tree$tip.label, drop = FALSE]
  depths <- diag(shared)
  patristic <- outer(depths, depths, "+") - 2 * shared
  patristic[patristic < 0 & patristic > -sqrt(.Machine$double.eps)] <- 0
  covariance <- (sigma2 / (2 * alpha)) * exp(-alpha * patristic) * (-expm1(-2 * alpha * shared))
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(tree$tip.label, tree$tip.label)
  if (any(!is.finite(covariance))) stop("OU covariance contains non-finite values.", call. = FALSE)
  covariance
}

fit_aic <- function(fit) {
  value <- fit$opt$aic
  if (is.null(value) || !is.finite(value)) stop("Model fit did not return a finite AIC.", call. = FALSE)
  value
}

select_model <- function(fits, force_model = NULL) {
  # PHYLOPHERE: explicit --perm_strategy override; NULL/auto keeps the AIC rule.
  if (!is.null(force_model) && nzchar(force_model)) {
    fm <- toupper(force_model)
    if (fm %in% c("BM", "OU")) return(fm)
  }
  bm <- fit_aic(fits$BM)
  ou <- fit_aic(fits$OU)
  if (ou + 2 < bm) "OU" else "BM"
}

analytical_s <- function(observed_difference, pair_variance) {
  tolerance <- 100 * .Machine$double.eps * max(abs(pair_variance), 1)
  pair_variance[pair_variance < 0 & pair_variance >= -tolerance] <- 0
  if (any(!is.finite(pair_variance) | pair_variance <= 0)) stop("Non-positive pairwise model variance.", call. = FALSE)
  z <- observed_difference / sqrt(pair_variance)
  -expm1(log(2) + stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
}

calculate_pairwise_scores <- function(values, tree, covariance_bm,
                                      covariance_ou, selected_model) {
  species <- tree$tip.label
  pairs <- utils::combn(seq_along(species), 2L)
  i <- pairs[1L, ]
  j <- pairs[2L, ]
  species_1 <- species[i]
  species_2 <- species[j]
  value_1 <- unname(values[species_1])
  value_2 <- unname(values[species_2])
  trait_difference <- abs(value_1 - value_2)
  patristic_matrix <- ape::cophenetic.phylo(tree)
  patristic_distance <- unname(patristic_matrix[cbind(species_1, species_2)])
  if (any(!is.finite(patristic_distance) | patristic_distance <= 0)) stop("Patristic distances must be positive and finite.", call. = FALSE)
  pair_variance <- function(covariance) {
    diag(covariance)[i] + diag(covariance)[j] - 2 * covariance[cbind(i, j)]
  }
  s_bm <- analytical_s(trait_difference, pair_variance(covariance_bm))
  s_ou <- analytical_s(trait_difference, pair_variance(covariance_ou))
  normalized_difference <- trait_difference / max(trait_difference)
  normalized_patristic <- patristic_distance / max(patristic_distance)
  selected_s <- if (selected_model == "BM") s_bm else s_ou
  final_score <- selected_s * normalized_difference / normalized_patristic
  out <- data.frame(
    Species1 = species_1, Species2 = species_2,
    TraitValue1 = value_1, TraitValue2 = value_2,
    PatristicDistance = patristic_distance,
    NormalizedTraitDifference = normalized_difference,
    NormalizedPatristicDistance = normalized_patristic,
    Sbm = s_bm, Sou = s_ou, FinalScore = final_score,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$FinalScore, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

model_fit_table <- function(fits, trait, selected_model) {
  data.frame(
    Trait = trait, Model = c("BM", "OU"),
    sigma2 = c(fits$BM$opt$sigsq, fits$OU$opt$sigsq),
    alpha = c(NA_real_, fits$OU$opt$alpha),
    AIC = c(fit_aic(fits$BM), fit_aic(fits$OU)),
    Selected = c("BM", "OU") == selected_model,
    stringsAsFactors = FALSE
  )
}

print.phylogenetic_shift_score_result <- function(x, ...) {
  cat("Analytical phylogenetic shift score\n")
  cat("Trait:", x$trait, "\n")
  cat("Species:", length(x$tree$tip.label), "\n")
  cat("Pairs:", nrow(x$scores), "\n")
  cat("Selected model:", x$selected_model, "\n")
  cat("Final-score range:", range(x$scores$FinalScore), "\n")
  invisible(x)
}
