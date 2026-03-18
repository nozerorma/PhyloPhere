# Plotting functions

# Resolve palette values at call time to avoid load-order issues
get_palette_values <- function(palette_name = NULL) {
  if (is.null(palette_name) && exists("color_palette", inherits = TRUE)) {
    palette_name <- get("color_palette", inherits = TRUE)
  }
  if (!is.character(palette_name) || !nzchar(palette_name)) {
    return(NULL)
  }
  if (!exists(palette_name, inherits = TRUE)) {
    return(NULL)
  }
  get(palette_name, inherits = TRUE)
}

resolve_taxa_palette <- function(taxa_values, palette_name = NULL, fallback_name = "fallback_palette") {
  taxa_values <- unique(as.character(stats::na.omit(taxa_values)))
  if (length(taxa_values) == 0) {
    return(NULL)
  }

  base_palette <- get_palette_values(palette_name)
  fallback_palette <- get_palette_values(fallback_name)
  if (is.null(fallback_palette) || length(fallback_palette) == 0) {
    fallback_palette <- grDevices::hcl.colors(max(20, length(taxa_values)), "Dynamic")
    names(fallback_palette) <- paste0("fallback_", seq_along(fallback_palette))
  }

  named_colors <- setNames(rep(NA_character_, length(taxa_values)), taxa_values)
  if (!is.null(base_palette) && length(base_palette) > 0) {
    matching_taxa <- intersect(taxa_values, names(base_palette))
    named_colors[matching_taxa] <- unname(base_palette[matching_taxa])
  }

  missing_taxa <- names(named_colors)[is.na(named_colors)]
  if (length(missing_taxa) > 0) {
    fallback_values <- unname(fallback_palette)
    fallback_idx <- vapply(
      missing_taxa,
      function(taxon_name) {
        ((sum(utf8ToInt(taxon_name)) - 1L) %% length(fallback_values)) + 1L
      },
      integer(1)
    )
    named_colors[missing_taxa] <- fallback_values[fallback_idx]
  }

  named_colors
}

clamp_value <- function(x, lower, upper) {
  max(lower, min(upper, x))
}

compute_ring_axis_limits <- function(values, min_limit = 1, headroom = 1.08) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    upper <- min_limit
  } else {
    upper <- max(values, na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) {
      upper <- min_limit
    } else {
      upper <- max(min_limit, upper * headroom)
    }
  }

  c(0, unname(upper))
}

species_plot_profile <- function(n_species, n_taxa = n_species, n_rings = 0L) {
  n_species <- max(1, as.numeric(n_species))
  n_taxa <- max(1, as.numeric(n_taxa))
  n_rings <- max(0, as.numeric(n_rings))

  contrast_height <- clamp_value(6 + 0.22 * n_taxa + 0.02 * n_species, 8, 22)
  contrast_width <- clamp_value(12 + 0.03 * n_species, 12, 18)
  violin_height <- clamp_value(6 + 0.24 * n_taxa + 0.015 * n_species, 8, 22)
  violin_width <- clamp_value(7 + 0.02 * n_species, 7, 12)

  tree_height <- clamp_value(10 + 0.11 * n_species, 12, 28)
  tree_width <- clamp_value(tree_height + 0.8 * n_rings + 2.5, 14, 34)
  diagnostic_size <- clamp_value(11 + 0.09 * n_species, 12, 28)

  list(
    n_species = n_species,
    n_taxa = n_taxa,
    n_rings = n_rings,
    contrast = list(
      width = contrast_width,
      height = contrast_height,
      point_size = clamp_value(5.0 - 0.018 * n_species, 2.0, 5.0),
      label_size = clamp_value(5.0 - 0.020 * n_species, 2.3, 5.0),
      axis_text_y = clamp_value(16.0 - 0.10 * n_taxa, 8.5, 16.0),
      title_size = clamp_value(20.0 - 0.08 * n_taxa, 14.0, 20.0),
      subtitle_size = clamp_value(12.0 - 0.04 * n_taxa, 9.0, 12.0),
      axis_title = clamp_value(15.0 - 0.04 * n_taxa, 11.0, 15.0),
      caption_size = clamp_value(12.0 - 0.03 * n_taxa, 9.0, 12.0),
      segment_size = clamp_value(0.32 - 0.002 * n_species, 0.12, 0.32),
      nudge_y = clamp_value(0.55 - 0.004 * n_species, 0.12, 0.55),
      nudge_x = clamp_value(0.010 - 0.00005 * n_species, 0.003, 0.010),
      force = clamp_value(1.0 - 0.004 * n_species, 0.25, 1.0),
      max_overlaps = clamp_value(round(40 - 0.25 * n_species), 8, 40),
      label_padding = grid::unit(clamp_value(0.22 - 0.0013 * n_species, 0.05, 0.22), "lines"),
      min_segment_length = clamp_value(0.06 - 0.0004 * n_species, 0.01, 0.06)
    ),
    violin = list(
      width = violin_width,
      height = violin_height,
      point_size = clamp_value(3.0 - 0.010 * n_species, 1.2, 3.0),
      stroke = clamp_value(0.8 - 0.003 * n_species, 0.25, 0.8),
      jitter_height = clamp_value(0.25 - 0.0015 * n_species, 0.06, 0.25),
      axis_text_y = clamp_value(17.0 - 0.11 * n_taxa, 8.5, 17.0),
      axis_title = clamp_value(17.0 - 0.05 * n_taxa, 11.0, 17.0)
    ),
    asr = list(
      width = clamp_value(14 + 0.08 * n_species, 14, 24),
      height = clamp_value(12 + 0.10 * n_species, 12, 28),
      fsize = clamp_value(1.8 - 0.011 * n_species, 0.55, 1.8),
      line_width = clamp_value(6.0 - 0.030 * n_species, 2.0, 6.0),
      main_cex = clamp_value(2.5 - 0.012 * n_species, 1.1, 2.5),
      colorbar_lwd = clamp_value(10.0 - 0.045 * n_species, 4.5, 10.0),
      colorbar_fsize = clamp_value(1.5 - 0.007 * n_species, 0.8, 1.5),
      tick_cex = clamp_value(1.5 - 0.007 * n_species, 0.7, 1.5),
      root_cex = clamp_value(1.5 - 0.007 * n_species, 0.7, 1.5),
      node_label_cex = clamp_value(1.5 - 0.008 * n_species, 0.65, 1.5),
      tip_symbol_cex = clamp_value(2.0 - 0.010 * n_species, 0.85, 2.0),
      tip_value_cex = clamp_value(1.5 - 0.007 * n_species, 0.7, 1.5),
      tip_symbol_offset = clamp_value(58.0 - 0.45 * n_species, 16.0, 58.0),
      tip_value_offset = clamp_value(40.0 - 0.30 * n_species, 10.0, 40.0),
      cumulative_offset = clamp_value(62.0 - 0.48 * n_species, 18.0, 62.0),
      segment_length = clamp_value(70.0 - 0.55 * n_species, 14.0, 70.0),
      segment_y = clamp_value(0.40 - 0.002 * n_species, 0.12, 0.40),
      legend_left = clamp_value(0.20 + 0.001 * n_species, 0.18, 0.30),
      xlim_right = clamp_value(2.00 + 0.010 * n_species, 2.00, 2.80),
      note_cex = clamp_value(1.2 - 0.004 * n_species, 0.75, 1.2)
    ),
    tree = list(
      width = tree_width,
      height = tree_height,
      diagnostic_size = diagnostic_size,
      tree_line_width = clamp_value(2.0 - 0.010 * n_species, 0.45, 2.0),
      node_text_size = clamp_value(6.0 - 0.030 * n_species, 1.8, 6.0),
      fruit_pwidth = clamp_value(0.75 - 0.004 * n_species, 0.24, 0.75),
      fruit_gap = clamp_value(0.10 + 0.0010 * n_species + 0.012 * pmax(0, n_rings - 2), 0.10, 0.22),
      secondary_offset = clamp_value(0.10 + 0.0010 * n_species + 0.012 * pmax(0, n_rings - 2), 0.10, 0.22),
      taxa_bar_pwidth = clamp_value(0.10 - 0.0003 * n_species, 0.05, 0.10),
      taxa_bar_offset = clamp_value(0.14 + 0.0014 * n_species + 0.012 * pmax(0, n_rings - 2), 0.14, 0.26),
      axis_text_size = clamp_value(7.2 - 0.020 * n_species, 3.6, 7.2),
      axis_nbreak = if (n_species > 80) 1 else 2,
      branch_label_size = clamp_value(8.5 - 0.030 * n_species, 4.5, 8.5),
      image_size = clamp_value(0.05 - 0.00025 * n_species, 0.02, 0.05),
      n_label_size = clamp_value(6.0 - 0.020 * n_species, 2.3, 6.0),
      n_label_nudge = clamp_value(8.4 + 0.60 * pmax(0, n_rings - 2) - 0.040 * n_species, 3.2, 9.8),
      asterisk_size = clamp_value(7.0 - 0.025 * n_species, 2.5, 7.0),
      asterisk_nudge = clamp_value(7.2 + 0.55 * pmax(0, n_rings - 2) - 0.032 * n_species, 2.8, 8.6),
      asterisk_y = clamp_value(0.18 - 0.0010 * n_species, 0.07, 0.18),
      cladelab_offset_phylopic = clamp_value(9.5 + 0.55 * n_rings + 0.03 * n_species, 9.5, 18.0),
      cladelab_offset_text = clamp_value(7.0 + 0.55 * n_rings + 0.025 * n_species, 7.0, 15.0),
      legend_title_size = clamp_value(15.0 - 0.030 * n_species, 10.0, 15.0),
      legend_text_size = clamp_value(13.0 - 0.025 * n_species, 8.0, 13.0),
      caption_size = clamp_value(15.0 - 0.025 * n_species, 9.0, 15.0)
    )
  )
}
