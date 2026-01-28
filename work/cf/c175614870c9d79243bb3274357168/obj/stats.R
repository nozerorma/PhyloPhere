
# Load required library
library(stats)

if (!exists("debug_log", inherits = TRUE)) {
  debug_log <- function(...) {
    msg <- sprintf(...)
    cat("[DEBUG] ", msg, "\n", sep = "")
  }
}

# Function to extract several statistical insights from the data
stats.f <- function(df) {
  p_trait_name <- if (exists("p_trait", inherits = TRUE)) p_trait else ""
  n_trait_name <- if (exists("n_trait", inherits = TRUE)) n_trait else ""
  tax_id <- if (exists("tax_id", inherits = TRUE)) tax_id else ""
  secondary_trait <- if (exists("secondary_trait", inherits = TRUE)) secondary_trait else ""
  branch_trait <- if (exists("branch_trait", inherits = TRUE)) branch_trait else ""
  taxa_col <- taxon_of_interest
  trait_col <- trait
  debug_log("stats.f rows = %d, trait = %s, taxon = %s, n_trait = %s, p_trait = %s, tax_id = %s", 
            nrow(df), trait_col, taxa_col, 
            ifelse(nzchar(n_trait_name), n_trait_name, "<none>"), 
            ifelse(nzchar(p_trait_name), p_trait_name, "<none>"),
            ifelse(nzchar(tax_id), tax_id, "<none>"),
            ifelse(nzchar(secondary_trait), secondary_trait, "<none>"),
            ifelse(nzchar(branch_trait), branch_trait, "<none>"))

  base_cols <- c("species", taxon_of_interest, trait)
  if (nzchar(n_trait_name) && n_trait_name %in% names(df)) {
    base_cols <- c(base_cols, n_trait_name)
  }
  if (nzchar(p_trait_name) && p_trait_name %in% names(df)) {
    base_cols <- c(base_cols, p_trait_name)
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
  g_mean <- mean(df_num[[trait]], na.rm = TRUE)
  g_median <- median(df_num[[trait]], na.rm = TRUE)
  g_sd <- sd(df_num[[trait]], na.rm = TRUE)
  g_q25 <- quantile(df_num[[trait]], 0.25, na.rm = TRUE)
  g_q75 <- quantile(df_num[[trait]], 0.75, na.rm = TRUE)
  g_iqr <- IQR(df_num[[trait]], na.rm = TRUE)

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
      g_q25 = g_q25,
      g_q75 = g_q75,
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
        (.data[[trait_col]] < g_q25) & (.data[[trait_col]] < g_median) ~ "low_extreme",
        (.data[[trait_col]] > g_q75) & (.data[[trait_col]] > g_median) ~ "high_extreme",
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
        (.data[[trait_col]] < taxa_q25) & (.data[[trait_col]] < taxa_median) ~ "low_extreme",
        (.data[[trait_col]] > taxa_q75) & (.data[[trait_col]] > taxa_median) ~ "high_extreme",
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
