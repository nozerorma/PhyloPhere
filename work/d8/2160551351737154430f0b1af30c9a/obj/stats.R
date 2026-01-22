
# Load required library
library(stats)

# Function to extract several statistical insights from the data
stats.f <- function(df) {
  taxa.sym <- rlang::sym(taxon_of_interest)
  trait.sym <- rlang::sym(trait)

  df <- df %>%
    dplyr::select(species, !!taxa.sym, !!trait.sym) %>%
    dplyr::filter(!is.na(!!trait.sym))

  if (nrow(df) < 4) {
    message("Skipping ", trait, " (fewer than 4 rows after NA filtering)")
    return(NULL)
  }

  if (!is_numeric_trait(df[[trait]])) {
    message("Skipping non-numeric trait: ", trait)
    return(NULL)
  }

  # Normalize to numeric once
  df_num <- df %>% dplyr::mutate(!!trait := as.numeric(!!trait.sym))

  # Global stats
  g_mean <- mean(df_num[[trait]], na.rm = TRUE)
  g_median <- median(df_num[[trait]], na.rm = TRUE)
  g_sd <- sd(df_num[[trait]], na.rm = TRUE)
  g_q25 <- quantile(df_num[[trait]], 0.25, na.rm = TRUE)
  g_q75 <- quantile(df_num[[trait]], 0.75, na.rm = TRUE)
  g_iqr <- IQR(df_num[[trait]], na.rm = TRUE)

  # Taxon stats separately, then join
  taxon_stats <- df_num %>%
    dplyr::group_by(!!taxa.sym) %>%
    dplyr::summarise(
      taxa_mean = mean(!!trait.sym, na.rm = TRUE),
      taxa_median = median(!!trait.sym, na.rm = TRUE),
      taxa_sd = sd(!!trait.sym, na.rm = TRUE),
      taxa_q25 = quantile(!!trait.sym, 0.25, na.rm = TRUE),
      taxa_q75 = quantile(!!trait.sym, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  df_num <- df_num %>%
    dplyr::left_join(taxon_stats, by = taxon_of_interest) %>%
    dplyr::mutate(
      value = !!trait.sym,
      g_mean = g_mean,
      g_median = g_median,
      g_sd = g_sd,
      g_q25 = g_q25,
      g_q75 = g_q75,
      outlier = dplyr::case_when(
        !!trait.sym < (g_q25 - 1.5 * g_iqr) ~ "low_outlier",
        !!trait.sym > (g_q75 + 1.5 * g_iqr) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      extreme_outlier = dplyr::case_when(
        !!trait.sym < (g_q25 - 3 * g_iqr) ~ "low_outlier",
        !!trait.sym > (g_q75 + 3 * g_iqr) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      global_label = dplyr::case_when(
        (!!trait.sym < g_q25) & (!!trait.sym < g_median) ~ "low_extreme",
        (!!trait.sym > g_q75) & (!!trait.sym > g_median) ~ "high_extreme",
        TRUE ~ "normal"
      )
    ) %>%
    dplyr::group_by(!!taxa.sym) %>%
    dplyr::mutate(
      taxa_outlier = dplyr::case_when(
        !!trait.sym < (taxa_q25 - 1.5 * IQR(!!trait.sym, na.rm = TRUE)) ~ "low_outlier",
        !!trait.sym > (taxa_q75 + 1.5 * IQR(!!trait.sym, na.rm = TRUE)) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      extreme_taxa_outlier = dplyr::case_when(
        !!trait.sym < (taxa_q25 - 3 * IQR(!!trait.sym, na.rm = TRUE)) ~ "low_outlier",
        !!trait.sym > (taxa_q75 + 3 * IQR(!!trait.sym, na.rm = TRUE)) ~ "high_outlier",
        TRUE ~ "normal"
      ),
      taxa_label = dplyr::case_when(
        (!!trait.sym < taxa_q25) & (!!trait.sym < taxa_median) ~ "low_extreme",
        (!!trait.sym > taxa_q75) & (!!trait.sym > taxa_median) ~ "high_extreme",
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