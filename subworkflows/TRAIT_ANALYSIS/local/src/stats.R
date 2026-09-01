
# Load required library
library(stats)

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Decide whether a trait vector is an ordinal fg/bg code rather than a
# continuous measurement. `trait_type` (from params / commons.R) forces it:
#   "ordinal"    -> always treat as coded
#   "continuous" -> never
#   ""/"auto"    -> infer: 2-5 distinct, all-integer levels => coded
# A genuine continuous trait (prevalence, body mass, ...) has many distinct
# non-integer values over the sampled species, so "auto" leaves it alone.
is_ordinal_trait <- function(trait_values, trait_type = "auto") {
  tt <- as.character(trait_type)
  if (length(tt) == 0 || is.na(tt[1]) || !nzchar(tt[1])) tt <- "auto"
  tt <- tolower(tt[1])
  if (tt == "ordinal") return(TRUE)
  if (tt == "continuous") return(FALSE)
  v <- trait_values[!is.na(trait_values)]
  u <- unique(v)
  length(u) >= 2 && length(u) <= 5 && all(u == round(u))
}

compute_trait_thresholds <- function(trait_values,
                                     discrete_method = if (exists("discrete_method", inherits = TRUE)) discrete_method else "decile",
                                     top_quantile = if (exists("top_quantile", inherits = TRUE)) as.numeric(top_quantile) else 0.90,
                                     bottom_quantile = if (exists("bottom_quantile", inherits = TRUE)) as.numeric(bottom_quantile) else 0.10) {
  trait_values <- as.numeric(trait_values)
  trait_values <- trait_values[!is.na(trait_values)]

  if (length(trait_values) == 0) {
    stop("compute_trait_thresholds() received no non-missing trait values.")
  }

  # trait_type: set by commons.R from params ("" / auto / ordinal / continuous).
  .trait_type <- if (exists("trait_type", inherits = TRUE)) get("trait_type", inherits = TRUE) else "auto"

  # ── Ordinal fg/bg code: the extremes ARE the top/bottom levels, no quantiles.
  # `quantile([mostly 0s], 0.9)` collapses to 0 and labels every species "top";
  # for a coded trait the highest level is foreground, the lowest is background,
  # any middle level is intermediate (kept in the tree, excluded from contrasts).
  if (is_ordinal_trait(trait_values, .trait_type)) {
    levels_sorted <- sort(unique(trait_values))
    lo <- levels_sorted[1]
    hi <- levels_sorted[length(levels_sorted)]
    return(list(
      discrete_method = "ordinal", top_quantile = NA_real_, bottom_quantile = NA_real_,
      trait_type = "ordinal", ordinal_levels = levels_sorted,
      g_mean = mean(trait_values), g_median = median(trait_values), g_sd = sd(trait_values),
      g_q10 = lo, g_q20 = lo, g_q25 = lo, g_q75 = hi, g_q80 = hi, g_q90 = hi,
      g_iqr = hi - lo,
      lower_thresh = unname(lo), upper_thresh = unname(hi)
    ))
  }

  g_mean <- mean(trait_values, na.rm = TRUE)
  g_median <- median(trait_values, na.rm = TRUE)
  g_sd <- sd(trait_values, na.rm = TRUE)
  g_q10 <- unname(quantile(trait_values, 0.10, na.rm = TRUE))
  g_q20 <- unname(quantile(trait_values, 0.20, na.rm = TRUE))
  g_q25 <- unname(quantile(trait_values, 0.25, na.rm = TRUE))
  g_q75 <- unname(quantile(trait_values, 0.75, na.rm = TRUE))
  g_q80 <- unname(quantile(trait_values, 0.80, na.rm = TRUE))
  g_q90 <- unname(quantile(trait_values, 0.90, na.rm = TRUE))
  g_iqr <- IQR(trait_values, na.rm = TRUE)

  lower_thresh <- switch(discrete_method,
    "quartile"      = g_q25,
    "quintile"      = g_q20,
    "decile"        = g_q10,
    "median_sd"     = g_median - g_sd,
    "parameterized" = unname(quantile(trait_values, bottom_quantile, na.rm = TRUE)),
    stop(sprintf(
      "Unknown discrete_method: '%s'. Choose one of: quartile, quintile, decile, median_sd, parameterized",
      discrete_method
    ))
  )
  upper_thresh <- switch(discrete_method,
    "quartile"      = g_q75,
    "quintile"      = g_q80,
    "decile"        = g_q90,
    "median_sd"     = g_median + g_sd,
    "parameterized" = unname(quantile(trait_values, top_quantile, na.rm = TRUE)),
    stop(sprintf(
      "Unknown discrete_method: '%s'. Choose one of: quartile, quintile, decile, median_sd, parameterized",
      discrete_method
    ))
  )

  list(
    discrete_method = discrete_method,
    top_quantile = top_quantile,
    bottom_quantile = bottom_quantile,
    g_mean = g_mean,
    g_median = g_median,
    g_sd = g_sd,
    g_q10 = g_q10,
    g_q20 = g_q20,
    g_q25 = g_q25,
    g_q75 = g_q75,
    g_q80 = g_q80,
    g_q90 = g_q90,
    g_iqr = g_iqr,
    lower_thresh = unname(lower_thresh),
    upper_thresh = unname(upper_thresh)
  )
}

# Function to extract several statistical insights from the data
stats.f <- function(df) {
  c_trait_name <- if (exists("c_trait", inherits = TRUE)) c_trait else ""
  n_trait_name <- if (exists("n_trait", inherits = TRUE)) n_trait else ""
  tax_id <- if (exists("tax_id", inherits = TRUE)) tax_id else ""
  secondary_trait <- if (exists("secondary_trait", inherits = TRUE)) secondary_trait else ""
  branch_trait <- if (exists("branch_trait", inherits = TRUE)) branch_trait else ""
  taxa_col <- taxon_of_interest
  trait_col <- trait
  debug_log("stats.f rows = %d, trait = %s, taxon = %s, n_trait = %s, c_trait = %s, tax_id = %s", 
            nrow(df), trait_col, taxa_col, 
            ifelse(nzchar(n_trait_name), n_trait_name, "<none>"), 
            ifelse(nzchar(c_trait_name), c_trait_name, "<none>"),
            ifelse(nzchar(tax_id), tax_id, "<none>"),
            ifelse(nzchar(secondary_trait), secondary_trait, "<none>"),
            ifelse(nzchar(branch_trait), branch_trait, "<none>"))

  base_cols <- c("species", taxon_of_interest, trait)
  if (nzchar(n_trait_name) && n_trait_name %in% names(df)) {
    base_cols <- c(base_cols, n_trait_name)
  }
  if (nzchar(c_trait_name) && c_trait_name %in% names(df)) {
    base_cols <- c(base_cols, c_trait_name)
  }
  if (nzchar(tax_id) && tax_id %in% names(df)) {
    base_cols <- c(base_cols, tax_id)
  }
  if (nzchar(secondary_trait) && secondary_trait %in% names(df)) {
    base_cols <- c(base_cols, secondary_trait)
  }
  if (nzchar(branch_trait) && branch_trait %in% names(df)) {
    base_cols <- c(base_cols, branch_trait)
  }

  debug_log("stats.f using columns: %s", paste(base_cols, collapse = ", "))

  df <- df %>%
    dplyr::select(dplyr::all_of(base_cols)) %>%
    dplyr::filter(!is.na(.data[[trait_col]]))
  debug_log("stats.f after NA filter rows = %d", nrow(df))

  if (nrow(df) < 4) {
    message("Skipping ", trait, " (fewer than 4 rows after NA filtering)")
    return(NULL)
  }

  # Normalize to numeric once
  df_num <- df %>% dplyr::mutate("{trait_col}" := as.numeric(.data[[trait_col]]))

  # Global stats
  disc_method <- if (exists("discrete_method", inherits = TRUE)) discrete_method else "decile"
  tq <- if (exists("top_quantile",    inherits = TRUE)) as.numeric(top_quantile)    else 0.90
  bq <- if (exists("bottom_quantile", inherits = TRUE)) as.numeric(bottom_quantile) else 0.10
  threshold_info <- compute_trait_thresholds(
    trait_values = df_num[[trait]],
    discrete_method = disc_method,
    top_quantile = tq,
    bottom_quantile = bq
  )
  g_mean <- threshold_info$g_mean
  g_median <- threshold_info$g_median
  g_sd <- threshold_info$g_sd
  g_q10 <- threshold_info$g_q10
  g_q20 <- threshold_info$g_q20
  g_q25 <- threshold_info$g_q25
  g_q75 <- threshold_info$g_q75
  g_q80 <- threshold_info$g_q80
  g_q90 <- threshold_info$g_q90
  g_iqr <- threshold_info$g_iqr
  lower_thresh <- threshold_info$lower_thresh
  upper_thresh <- threshold_info$upper_thresh
  debug_log("stats.f: disc_method = %s, lower_thresh = %.4f, upper_thresh = %.4f",
            disc_method, lower_thresh, upper_thresh)

  # Taxon stats separately, then join
  taxon_stats <- df_num %>%
    dplyr::group_by(.data[[taxa_col]]) %>%
    dplyr::summarise(
      taxa_mean = mean(.data[[trait_col]], na.rm = TRUE),
      taxa_median = median(.data[[trait_col]], na.rm = TRUE),
      taxa_sd = sd(.data[[trait_col]], na.rm = TRUE),
      taxa_q25 = quantile(.data[[trait_col]], 0.25, na.rm = TRUE),
      taxa_q75 = quantile(.data[[trait_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  df_num <- df_num %>%
    dplyr::left_join(taxon_stats, by = taxon_of_interest) %>%
    dplyr::mutate(
      value = .data[[trait_col]],
      g_mean = g_mean,
      g_median = g_median,
      g_sd = g_sd,
      g_q10 = g_q10,
      g_q20 = g_q20,
      g_q25 = g_q25,
      g_q75 = g_q75,
      g_q80 = g_q80,
      g_q90 = g_q90,
      lower_thresh = lower_thresh,
      upper_thresh = upper_thresh,
      outlier = dplyr::case_when(
        .data[[trait_col]] < (g_q25 - 1.5 * g_iqr) ~ "low_outlier",
        .data[[trait_col]] > (g_q75 + 1.5 * g_iqr) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      extreme_outlier = dplyr::case_when(
        .data[[trait_col]] < (g_q25 - 3 * g_iqr) ~ "low_outlier",
        .data[[trait_col]] > (g_q75 + 3 * g_iqr) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      global_label = dplyr::case_when(
        # `<= g_median` / `>= g_median` (was `<` / `>`): the strict form makes a
        # minority-foreground binary/bimodal trait un-labelable — median == the
        # lower threshold == 0, so no species can be `< g_median`, and FADE's
        # EXTRACT_EXTREME_SPECIES then finds no low_extreme set. No-op for a
        # continuous trait (upper_thresh >= median >= lower_thresh always).
        # Mirrors the same relaxation in lean_contrast_selector.R.
        (.data[[trait_col]] <= lower_thresh) & (.data[[trait_col]] <= g_median) ~ "low_extreme",
        (.data[[trait_col]] >= upper_thresh) & (.data[[trait_col]] >= g_median) ~ "high_extreme",
        TRUE ~ "normal"
      )
    ) %>%
    dplyr::group_by(.data[[taxa_col]]) %>%
    dplyr::mutate(
      taxa_outlier = dplyr::case_when(
        .data[[trait_col]] < (taxa_q25 - 1.5 * IQR(.data[[trait_col]], na.rm = TRUE)) ~ "low_outlier",
        .data[[trait_col]] > (taxa_q75 + 1.5 * IQR(.data[[trait_col]], na.rm = TRUE)) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      extreme_taxa_outlier = dplyr::case_when(
        .data[[trait_col]] < (taxa_q25 - 3 * IQR(.data[[trait_col]], na.rm = TRUE)) ~ "low_outlier",
        .data[[trait_col]] > (taxa_q75 + 3 * IQR(.data[[trait_col]], na.rm = TRUE)) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      taxa_label = dplyr::case_when(
        # median-guard relaxed to `<=` / `>=` as in global_label above
        (.data[[trait_col]] < taxa_q25) & (.data[[trait_col]] <= taxa_median) ~ "low_extreme",
        (.data[[trait_col]] > taxa_q75) & (.data[[trait_col]] >= taxa_median) ~ "high_extreme",
        TRUE ~ "normal"
      )
    ) %>%
    dplyr::ungroup()
  df_num
}

# Robust scaling
robust_scale <- function(x) {
  (x - median(x)) / IQR(x)
}

# Find outliers
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5 * IQR(x) | x > quantile(x, .75) + 1.5 * IQR(x))
}

# ----------------------------------------
# Confidence Interval Functions
# ----------------------------------------

#' Check if two confidence intervals overlap
#' 
#' @param lb1 Lower bound of first interval
#' @param ub1 Upper bound of first interval
#' @param lb2 Lower bound of second interval
#' @param ub2 Upper bound of second interval
#' @return Logical indicating whether intervals overlap
ci_overlap <- function(lb1, ub1, lb2, ub2) {
  (lb1 <= ub2) & (lb2 <= ub1)
}
