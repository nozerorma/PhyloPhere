# Observed composition/selection contract. Functions take explicit inputs so the
# two reports share one prepared tree, fit, score table and candidate gate.
contrast_integer <- function(x, name, minimum = 0) {
  n <- suppressWarnings(as.numeric(x))
  if (length(n) != 1L || !is.finite(n) || n != floor(n) ||
      n < minimum || n > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, ".")
  as.integer(n)
}

validate_contrast_inputs <- function(df, trait, n_trait, c_trait, params) {
  for (column in c(n_trait, c_trait)) {
    if (nzchar(column) && !column %in% names(df))
      stop("Configured count column not found: ", column)
  }
  if (nzchar(c_trait) && !nzchar(n_trait)) stop("c_trait requires n_trait.")
  if (nzchar(n_trait) && nzchar(c_trait)) {
    if (length(unique(c(trait, n_trait, c_trait))) != 3L)
      stop("traitname, n_trait and c_trait must name distinct columns.")
    n <- df[[n_trait]]; k <- df[[c_trait]]
    bad <- !is.na(n) & !is.na(k) &
      (!is.finite(n) | !is.finite(k) | n <= 0 | k < 0 | k > n |
       n != floor(n) | k != floor(k))
    if (any(bad)) stop("Invalid integer counts (require n > 0 and 0 <= c <= n): ",
                       paste(df$species[bad], collapse = ", "))
  }
  mode <- if (is.null(params$trait_type) || !nzchar(params$trait_type)) "auto" else tolower(params$trait_type)
  if (!mode %in% c("auto", "ordinal", "continuous")) stop("Invalid trait_type: ", mode)
  strategy <- if (is.null(params$perm_strategy)) "best_model" else tolower(params$perm_strategy)
  if (!strategy %in% c("auto", "best_model", "bm", "ou")) stop("Invalid perm_strategy: ", strategy)
  pct <- suppressWarnings(as.numeric(params$pss_top_pct))
  if (length(pct) != 1L || !is.finite(pct) || pct < 0 || pct > 1)
    stop("pss_top_pct must be a finite fraction from 0 to 1 (0 disables the gate).")
  contrast_integer(params$seed, "seed")
  if (!is.null(params$max_contrasts)) contrast_integer(params$max_contrasts, "max_contrasts")
  if (!is.null(params$max_fop)) contrast_integer(params$max_fop, "max_fop", 1)
  invisible(TRUE)
}

contrast_provenance <- function(df, tree, trait, n_trait, c_trait,
                                top_pct, ordinal, center_gate, force_model, seed) {
  list(version = 1L, data = df[, unique(c("species", trait, n_trait[nzchar(n_trait)],
       c_trait[nzchar(c_trait)])), drop = FALSE], tree = tree,
       top_pct = as.numeric(top_pct), ordinal = ordinal, center_gate = center_gate,
       force_model = force_model, seed = as.integer(seed))
}

prepare_observed_contrasts <- function(ci, tree, top_pct, ordinal, center_gate,
                                       force_model = NULL) {
  v <- setNames(ci$trait, ci$species)
  use_ci <- any(!is.na(ci$n_population))
  bounds <- function(column) if (use_ci) setNames(ci[[column]], ci$species) else NULL
  lb <- bounds("trait_ci_lb"); ub <- bounds("trait_ci_ub"); n <- bounds("n_population")
  # BM/OU fitting needs >= 3 tips and variation. These are valid zero-pair
  # outcomes, unlike an optimizer/model failure on otherwise usable data.
  if (length(v) < 3L || length(unique(v)) < 2L) {
    D <- if (length(v) > 1L) ape::cophenetic.phylo(tree) else
      matrix(0, 1L, 1L, dimnames = list(names(v), names(v)))
    ctx <- list(tree = tree, trait_vec = v, D = D, selected_model = NA_character_)
    return(list(ctx = ctx, candidates = list(cand_df = NULL, mode = "empty",
      reason = "fewer than three usable species or no trait variation"), scores = NULL))
  }
  ctx <- tryCatch(selection_context(v, tree, force_model), error = function(e)
    stop("PSS model fitting failed: ", conditionMessage(e), call. = FALSE))
  sc <- calculate_pairwise_scores(ctx$trait_vec, ctx$tree, ctx$cov_bm,
                                 ctx$cov_ou, ctx$selected_model)
  if (use_ci) {
    higher <- ifelse(sc$TraitValue1 > sc$TraitValue2, sc$Species1, sc$Species2)
    lower <- ifelse(sc$TraitValue1 > sc$TraitValue2, sc$Species2, sc$Species1)
    reversed <- sc$TraitValue1 != sc$TraitValue2 & ub[higher] < lb[lower]
    if (any(reversed)) stop("Count intervals conflict with phenotype orientation: ",
      paste(paste(higher[reversed], lower[reversed], sep = "/"), collapse = ", "))
  }
  cc <- lean_candidate_df(ctx$trait_vec, ctx$D, 1L, ctx$tree, ctx$cov_bm,
    ctx$cov_ou, ctx$selected_model, lb, ub, top_pct, ordinal, n, center_gate, scores = sc)
  list(ctx = ctx, candidates = cc, scores = sc)
}

read_observed_contrasts <- function(path, provenance) {
  if (!file.exists(path)) stop("Missing observed selection inputs; rerun composition: ", path)
  bundle <- readRDS(path)
  if (!identical(bundle$provenance, provenance))
    stop("Observed selection provenance mismatch; rerun composition.")
  bundle
}

gen_traitfiles <- function(df) {
  # cluster is the domain id used by the scoring metadata.
  pair <- if ("cluster" %in% names(df)) as.integer(df$cluster) else seq_len(nrow(df))
  result <- rbind(data.frame(species = df$species1, contrast_group = rep.int(1L, nrow(df)), pair = pair),
                  data.frame(species = df$species2, contrast_group = rep.int(0L, nrow(df)), pair = pair))
  result[order(result$contrast_group, result$pair), , drop = FALSE]
}
