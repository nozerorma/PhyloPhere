# Plotting functions

# Resolve palette values at call time to avoid load-order issues
get_palette_values <- function() {
  if (!exists("color.sym", inherits = TRUE)) {
    return(NULL)
  }
  palette_name <- rlang::as_string(color.sym)
  if (!exists(palette_name, inherits = TRUE)) {
    return(NULL)
  }
  get(palette_name, inherits = TRUE)
}

## 1.Dataset_exploration
taxa_plot.f <- function(df){
  max_count <- max(df$count, na.rm = TRUE)
  palette_values <- get_palette_values()
  fill_scale <- if (is.null(palette_values)) {
    ggplot2::scale_fill_discrete()
  } else {
    ggplot2::scale_fill_manual(values = palette_values, breaks = NULL)
  }
  ggplot2::ggplot(df, 
                aes(x = !!taxa.sym, y = count, fill = !!taxa.sym)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), 
            vjust = 0.5, hjust = -1.5, 
            color = "black", fontface = "bold", 
            size = 6) +
  labs(x = paste(capitalized_taxon, "-", clade_name), 
       y = "Number of species", 
       title = paste(capitalized_taxon, "distribution in", clade_name)) +
  fill_scale + 
  scale_y_continuous(limits = c(0, max_count+2), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(hjust = 1, size = 18, face = "bold", color = "black"),
    axis.text.x = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20), face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, hjust = 0.5, margin = margin(t = 20), face = "bold", color = "black"),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  coord_flip()
}

## 2.Dataset_exploration
### Contrast plot function
contrast_plot.f <- function(df) {
  stopifnot(all(c("species", taxon_of_interest, trait) %in% names(df)))
  palette_values <- get_palette_values()
  fill_scale <- if (is.null(palette_values)) {
    ggplot2::scale_fill_discrete()
  } else {
    ggplot2::scale_fill_manual(values = palette_values)
  }
  color_scale <- if (is.null(palette_values)) {
    ggplot2::scale_color_discrete()
  } else {
    ggplot2::scale_color_manual(values = palette_values)
  }
  
  df <- df %>%
    dplyr::mutate(value = !!trait.sym, category = global_label, taxa = !!taxa.sym) %>%
    dplyr::mutate(common_name = gsub("_", " ", species))
  
  if(has.N == TRUE){
    df <- df %>%   
      mutate(N = !!N.sym)
  }
  
  # Ensure species order matches phylogeny
  ordered_species <- pruned_tree$tip.label
  df$Species <- factor(df$species, levels = ordered_species)
  
  # Define ordered taxa for plotting
  ordered_taxa <- unique(df$taxa)

  select_cols <- c(
    "taxa", "species", "common_name", "value", "category",
    "outlier", "extreme_outlier", "g_median", "taxa_median"
  )
  if (isTRUE(has.N)) {
    select_cols <- c(select_cols, "N")
  }

  df <- df %>%
    dplyr::arrange(taxa, .by_group = TRUE) %>%
    dplyr::select(dplyr::all_of(select_cols))
  
  global_median <- df$g_median * 100
  
  ggplot(df, aes(x = value * 100, y = taxa)) +
    ggplot2::scale_y_discrete(limits = ordered_taxa) +
    ggplot2::geom_point(
      data = subset(df, category == "normal"),
      aes(shape = category, color = taxa), size = 5
    ) +
    ggplot2::geom_point(
      data = subset(df, category %in% c("low_extreme", "high_extreme")),
      aes(shape = category, color = taxa, fill = taxa), size = 5
    ) +
    ggplot2::scale_shape_manual(values = c("low_extreme" = 15, "normal" = 1, "high_extreme" = 17)) +
    fill_scale +
    color_scale +
    ggplot2::geom_vline(xintercept = global_median, linetype = "longdash", linewidth = 1.2, color = "salmon3") +
    ggplot2::stat_summary(fun = median, geom = "errorbar",
                          aes(xmax = after_stat(x), xmin = after_stat(x), y = taxa),
                          linewidth = 1.2, color = "skyblue4", alpha = 0.8) +
    ggplot2::labs(x = paste(capitalized_trait, "(%)"), 
                  y = paste(capitalized_taxon, "-", clade_name), 
                  title = paste("Trait differential plot for trait:", capitalized_trait, "in", clade_name)) +
    ggrepel::geom_text_repel(
      data = subset(df, category %in% c("low_extreme", "high_extreme")),
      aes(label = species), size = 5, hjust = 0, vjust = 0, family = "Inter",
      segment.size = 0.3, nudge_x = 0.010, nudge_y = 0.5, max.overlaps = Inf,
      force = 1, direction = "y", max.iter = 10000, label.padding = 20,
      min.segment.length = 0.06, seed = 2000, show.legend = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(size = 15, hjust = 0.5, margin = ggplot2::margin(t = 20, b = 20)),
      axis.title.y = ggplot2::element_text(size = 15, angle = 90, hjust = 0.5, margin = ggplot2::margin(l = 20, r = 20)),
      legend.text = ggplot2::element_text(size = 15),
      legend.position = "right"
    )
}
