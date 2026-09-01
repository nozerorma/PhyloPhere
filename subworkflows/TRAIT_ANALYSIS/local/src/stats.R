
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

# Descriptive global summary of a trait vector: location, spread, and the
# quantiles used for IQR outlier fences and summary tables. This is pure
# description — the fg/bg (foreground/background) partition that drives contrast
# selection and FADE is produced downstream by 4.Independent_contrasts.Rmd from
# the phylogenetically-independent selected pairs (PSS for continuous traits,
# extreme coded levels for ordinal), never here.
compute_trait_summary <- function(trait_values) {
  trait_values <- as.numeric(trait_values)
  trait_values <- trait_values[!is.na(trait_values)]

  if (length(trait_values) == 0) {
    stop("compute_trait_summary() received no non-missing trait values.")
  }

  list(
    g_mean   = mean(trait_values),
    g_median = median(trait_values),
    g_sd     = sd(trait_values),
    g_q10    = unname(quantile(trait_values, 0.10)),
    g_q20    = unname(quantile(trait_values, 0.20)),
    g_q25    = unname(quantile(trait_values, 0.25)),
    g_q75    = unname(quantile(trait_values, 0.75)),
    g_q80    = unname(quantile(trait_values, 0.80)),
    g_q90    = unname(quantile(trait_values, 0.90)),
    g_iqr    = IQR(trait_values)
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

  # Global descriptive stats (location, spread, quantiles for outlier fences).
  # No fg/bg partition here — 4.Independent_contrasts.Rmd owns that.
  trait_summary <- compute_trait_summary(df_num[[trait]])
  g_mean <- trait_summary$g_mean
  g_median <- trait_summary$g_median
  g_sd <- trait_summary$g_sd
  g_q10 <- trait_summary$g_q10
  g_q20 <- trait_summary$g_q20
  g_q25 <- trait_summary$g_q25
  g_q75 <- trait_summary$g_q75
  g_q80 <- trait_summary$g_q80
  g_q90 <- trait_summary$g_q90
  g_iqr <- trait_summary$g_iqr
  debug_log("stats.f: g_median = %.4f, g_iqr = %.4f", g_median, g_iqr)

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
      outlier = dplyr::case_when(
        .data[[trait_col]] < (g_q25 - 1.5 * g_iqr) ~ "low_outlier",
        .data[[trait_col]] > (g_q75 + 1.5 * g_iqr) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      extreme_outlier = dplyr::case_when(
        .data[[trait_col]] < (g_q25 - 3 * g_iqr) ~ "low_outlier",
        .data[[trait_col]] > (g_q75 + 3 * g_iqr) ~ "high_outlier",
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
        # descriptive per-taxon tail flag for exploration plots only; the
        # median guard keeps a minority-tail trait labelable (see compute_trait_summary)
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
