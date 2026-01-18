# Function to extract several statistical insights from the data
stats.f <- function(df) {
  
  df <- df %>%
    select(species, !!taxa.sym, !!trait.sym) %>%
    filter(!is.na(!!trait.sym))
  
  if (nrow(df) < 4) {
    message("Skipping ", trait, " (fewer than 4 rows after NA filtering)")
    return(NULL)
  }
  
  # Numeric branch
  if (is_numeric_trait(df[[trait]])) {
    message("Numeric trait: ", trait)
    
    # Normalize to numeric once
    df_num <- df %>% mutate(!!trait := as.numeric(!!trait.sym))
    
    # Global stats
    g_mean   <- mean(df_num[[trait]],   na.rm = TRUE)
    g_median <- median(df_num[[trait]], na.rm = TRUE)
    g_sd     <- sd(df_num[[trait]],     na.rm = TRUE)
    
    # Extremes (defined by quantiles)
    g_q25 <- quantile(df_num[[trait]], 0.25, na.rm = TRUE)
    g_q75 <- quantile(df_num[[trait]], 0.75, na.rm = TRUE)
    
    # Taxon stats separately, then join
    taxon_stats <- df_num %>%
      group_by(!!taxa.sym) %>%
      summarise(
        taxa_mean   = mean(!!trait.sym, na.rm = TRUE),
        taxa_median = median(!!trait.sym, na.rm = TRUE),
        taxa_sd     = sd(!!trait.sym, na.rm = TRUE),
        .groups = "drop"
      )
    
    df_num <- df_num %>%
      left_join(taxon_stats, by = taxon_of_interest) %>%
      mutate(
        # Global outlier by 1.5*SD rule
        outlier = case_when(
          !!trait.sym > (g_mean - 1.5 * g_sd) ~ "low_outlier",
          !!trait.sym > (g_mean + 1.5 * g_sd) ~ "high_outlier",
          TRUE ~ "normal"
        ) %>%
          # Extreme global outlier by 3*SD rule
          extreme_outlier = case_when(
            !!trait.sym > (g_mean - 3 * g_sd) ~ "low_outlier",
            !!trait.sym > (g_mean + 3 * g_sd) ~ "high_outlier",
            TRUE ~ "normal"
          ) %>%        
          # Phenotypic labels using precomputed global quantiles/median
          global_label = case_when(
            (!!trait.sym < g_q25) & (!!trait.sym < g_median) ~ "low_extreme",
            (!!trait.sym > g_q75) & (!!trait.sym > g_median) ~ "high_extreme",
            TRUE ~ "normal"
          )
      ) %>%
      # family labels: compute per-family quantiles; do it within group
      group_by(!!taxa.sym) %>%
      mutate(
        taxa_q25   = quantile(!!trait.sym, 0.25, na.rm = TRUE),
        taxa_q75   = quantile(!!trait.sym, 0.75, na.rm = TRUE),
        # Global outlier by 1.5*SD rule
        taxa_outlier = case_when(
          !!trait.sym > (taxa_mean - 1.5 * taxa_sd) ~ "low_outlier",
          !!trait.sym > (taxa_mean + 1.5 * taxa_sd) ~ "high_outlier",
          TRUE ~ "normal"
        ) %>%
          # Extreme global outlier by 3*SD rule
          extreme_taxa_outlier = case_when(
            !!trait.sym > (taxa_mean - 3 * taxa_sd) ~ "low_outlier",
            !!trait.sym > (taxa_mean + 3 * taxa_sd) ~ "high_outlier",
            TRUE ~ "normal"
          ) %>%  
          taxa_label = case_when(
            (!!trait.sym < taxa_q25) & (!!trait.sym < taxa_median) ~ "low_extreme",
            (!!trait.sym > taxa_q75) & (!!trait.sym > taxa_median) ~ "high_extreme",
            TRUE ~ "normal"
          )
      ) %>%
      ungroup() %>%
      # mutate(value = !!t_sym) %>%
      select(-taxa_q25, -taxa_q75) # drop helper columns
    
    return(df_num)
  }
}
