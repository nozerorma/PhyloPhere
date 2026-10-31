# report_plots.R — Shared contrast, violin and ASR plots for reports 0, 2 and 4.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# Sourced locally by reports 0, 2 and 4 so labels and palettes use the render environment.

contrast_plot.f <- function(df, palette_values = NULL, point_shape = 16) {
  stopifnot(all(c("species", taxon_of_interest, trait) %in% names(df)))
  plot_profile <- species_plot_profile(
    n_species = dplyr::n_distinct(df$species),
    n_taxa = dplyr::n_distinct(df$taxa)
  )
  contrast_dims <- plot_profile$contrast
  if (is.null(palette_values)) {
    palette_values <- resolve_taxa_palette(df$taxa, palette_name = color_palette)
  }
  debug_log("contrast_plot.f rows = %d, trait = %s, taxon = %s", nrow(df), trait, taxon_of_interest)
  debug_log("contrast_plot.f columns: %s", paste(names(df), collapse = ", "))
  debug_log("color_palette = %s", color_palette)
  debug_log("palette_values is NULL: %s", is.null(palette_values))
  debug_log("palette_values length: %d", ifelse(is.null(palette_values), 0, length(palette_values)))
  if (!is.null(palette_values)) {
    debug_log("palette_values: %s", paste(palette_values, collapse = ", "))
  }
  
  global_median <- df$g_median[1]

  # One uniform point layer per species.
  p <- ggplot(df, aes(x = value, y = taxa)) +
    ggplot2::scale_y_discrete(limits = rev(ordered_taxa)) +
    ggplot2::geom_point(
      aes(color = taxa), shape = point_shape, size = contrast_dims$point_size
    )

  # breaks = 0 matches no level, which hides the legends and keeps the colors.
  if (!is.null(palette_values) && length(palette_values) > 0) {
    p <- p +
      ggplot2::scale_fill_manual(values = palette_values, breaks = 0) +
      ggplot2::scale_color_manual(values = palette_values, breaks = 0)
  }

  # Reference marks: the global median and the per-taxon medians.
  p <- p +
    ggplot2::geom_vline(xintercept = global_median, linetype = "longdash", linewidth = 1.2, color = tableau_10["vermillion"]) +
    ggplot2::stat_summary(fun = median, geom = "errorbar",
                          aes(xmax = after_stat(x), xmin = after_stat(x), y = taxa),
                          linewidth = 1.2, color = tableau_10["blue"], alpha = 0.8) +
    ggplot2::labs(x = capitalized_trait,
                  y = paste(capitalized_taxon, "-", capitalized_clade),
                  title = paste("Trait distribution for trait:", capitalized_trait, "in", clade_name),
                  caption = "Global median: longdash red line. Per-taxa median: blue error bars.") +
    report_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = contrast_dims$title_size, hjust = 0, margin = ggplot2::margin(t = 20, b = 6), face = "bold"),
      plot.subtitle = ggplot2::element_text(size = contrast_dims$subtitle_size, hjust = 0, color = "grey40", margin = ggplot2::margin(b = 10)),
      axis.text.y = ggplot2::element_text(size = contrast_dims$axis_text_y),
      axis.title.x = ggplot2::element_text(size = contrast_dims$axis_title, margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(size = contrast_dims$axis_title, angle = 90, margin = ggplot2::margin(l = 20, r = 20)),
      plot.caption = ggplot2::element_text(size = contrast_dims$caption_size, hjust = 0, margin = ggplot2::margin(t = 20)),
      legend.text = ggplot2::element_text(size = 15),
      legend.position = "bottom"
    )
  attr(p, "plot_dims") <- contrast_dims[c("width", "height")]
  p
}


violin_extremes.f <- function(df, trait, taxon_of_interest, point_shape = 16) {
  stopifnot(all(c("species", taxon_of_interest, trait) %in% names(df)))
  plot_profile <- species_plot_profile(
    n_species = dplyr::n_distinct(df$species),
    n_taxa = dplyr::n_distinct(df$taxa)
  )
  violin_dims <- plot_profile$violin
  palette_values <- resolve_taxa_palette(df$taxa, palette_name = color_palette)
  debug_log("violin_extremes.f rows = %d, trait = %s, taxon = %s", nrow(df), trait, taxon_of_interest)
  fill_scale <- if (is.null(palette_values)) {
    ggplot2::scale_fill_discrete()
  } else {
    ggplot2::scale_fill_manual(values = palette_values)
  }
  dark_palette_values <- resolve_taxa_palette(df$taxa, paste0("dark_", color_palette))
  fill_scale_dark <- if (is.null(dark_palette_values)) {
    fill_scale
  } else {
    ggplot2::scale_fill_manual(values = dark_palette_values)
  }
  color_scale_dark <- if (is.null(dark_palette_values)) {
    ggplot2::scale_color_discrete()
  } else {
    ggplot2::scale_color_manual(values = dark_palette_values)
  }

  # Species in tree order.
  ordered_species <- pruned_tree$tip.label
  df$species <- factor(df$species, levels = ordered_species)
  debug_log("violin_extremes.f ordered_species = %d", length(ordered_species))
  
  # Taxa in tree order: the global ordered_taxa when it exists, else rebuilt from the tree.
  if (!exists("ordered_taxa", inherits = TRUE) || length(ordered_taxa) == 0) {
    ordered_taxa <- unique(pruned_tree$tip.label %>%
      sapply(function(sp) df$taxa[df$species == sp]) %>%
      na.omit())
    if (length(ordered_taxa) == 0) {
      ordered_taxa <- unique(df$taxa[!is.na(df$taxa)])
    }
  }
  
  median_value <- median(df$value, na.rm = TRUE)
  
  plot <- ggplot2::ggplot(data = df, aes(x = value, y = taxa)) +
    ggplot2::scale_y_discrete(limits = rev(ordered_taxa)) +
    ggplot2::geom_violin(
      aes(fill = taxa), color = "black", width = 0.5, trim = TRUE, scale = "width", alpha = 0.15
    ) +
    fill_scale_dark +
    ggplot2::geom_point(
      data = subset(df, outlier == "normal"),
      aes(color = taxa), shape = point_shape,
      size = violin_dims$point_size, position = ggplot2::position_jitter(height = violin_dims$jitter_height)
    ) +
    ggplot2::geom_point(
      data = subset(df, outlier != "normal"),
      aes(color = taxa, fill = taxa),
      shape = 21, size = violin_dims$point_size, stroke = violin_dims$stroke,
      position = ggplot2::position_jitter(height = violin_dims$jitter_height)
    ) +
    color_scale_dark +
    ggplot2::geom_vline(
      aes(xintercept = median_value), linewidth = 1, color = tableau_10["vermillion"], alpha = 0.5
    ) +
    ggplot2::stat_summary(
      fun = median, geom = "errorbar",
      aes(xmax = after_stat(x), xmin = after_stat(x)),
      linewidth = 1.5, color = "black", alpha = 0.6, width = 0.75
    ) +
    ggplot2::labs(x = capitalized_trait) +
    report_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = violin_dims$axis_text_y, hjust = 0.5, face = "bold"),
      axis.title.x = ggplot2::element_text(size = violin_dims$axis_title, face = "bold", margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 15),
      legend.position = "none"
    )

  attr(plot, "plot_dims") <- violin_dims[c("width", "height")]
  plot
}


asr_tree.f <- function(df, tree, trait, taxon_of_interest,
                       fit_fun = geiger::fitContinuous, ace_fun = ape::ace) {
  
  # Trait vector, named by species, without missing or non-finite values.
  trait_vec <- df[[trait]]
  names(trait_vec) <- df$species
  trait_vec <- trait_vec[!is.na(trait_vec) & is.finite(trait_vec)]
  if(length(trait_vec) < 5) stop("ASR model comparison requires at least five species.")
  debug_log("asr_tree.f input rows = %d, trait_vec = %d", nrow(df), length(trait_vec))
  

  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, names(trait_vec)))
  # ape::drop.tip() can leave a $root.edge on the pruned tree even when the source Newick
  # had none (when the root keeps a single child, its branch length is folded in as a root
  # edge). The root edge is removed before the fit, rescaling and ASR steps.
  tree$root.edge <- NULL
  tree <- ape::reorder.phylo(tree, "cladewise")
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  asr_dims <- species_plot_profile(n_species = Ntip, n_taxa = dplyr::n_distinct(df[[taxon_of_interest]]))$asr
  debug_log("asr_tree.f tree tips after prune = %d, nodes = %d", Ntip, Nnode)
  node_indices <- (Ntip+1):(Ntip+Nnode)

  # Trait vector and data frame in the tip order of the pruned tree.
  trait_vec <- trait_vec[tree$tip.label]

  df <- df[df$species %in% tree$tip.label, ]
  df$species <- factor(df$species, levels = tree$tip.label)

  debug_log("asr_tree.f trait_vec after prune = %d", length(trait_vec))
  debug_log("asr_tree.f df after prune = %d", nrow(df))
  
  # Fit BM and OU; sigma2 is the rate and alpha the OU pull (0 under BM).
  fit_model <- function(tree, trait, model = c("BM", "OU")) {
    model <- match.arg(model)
    fit <- if (model == "OU") {
      # alpha upper bound scaled to tree height, as in pss_core.R::ou_alpha_bounds
      suppressWarnings(fit_fun(
        tree, trait, model = "OU", ncores = 1,
        bounds = list(alpha = c(exp(-500), 100 / max(ape::node.depth.edgelength(tree))))
      ))
    } else {
      suppressWarnings(fit_fun(tree, trait, model = "BM", ncores = 1))
    }
    aicc <- fit$opt$aicc
    list(
      lnL = fit$opt$lnL,
      sigma2 = fit$opt$sigsq,
      alpha = if (model == "OU") fit$opt$alpha else 0,
      aicc = aicc,
      raw = fit
    )
  }
  
  fit_bm <- fit_model(tree, trait_vec, "BM")
  fit_ou <- fit_model(tree, trait_vec, "OU")
  debug_log("asr_tree.f fit_bm aicc = %.3f, fit_ou aicc = %.3f", fit_bm$aicc, fit_ou$aicc)
  
  delta_aicc <- fit_bm$aicc - fit_ou$aicc
  chosen_model <- if (delta_aicc >= 2.0 && !is.null(fit_ou$alpha) && fit_ou$alpha > 1e-4) "OU" else "BM"
  debug_log("asr_tree.f chosen_model = %s", chosen_model)
  
  # Under OU the tree is rescaled with the fitted alpha before the ASR.
  if (chosen_model == "OU") {
    alpha_hat <- fit_ou$alpha
    debug_log("asr_tree.f alpha_hat = %.4f", alpha_hat)
    sigma2_used <- fit_ou$sigma2
    rescaled_tree <- phytools::rescale(tree, "OU", alpha = alpha_hat)
  } else {
    sigma2_used <- fit_bm$sigma2
    rescaled_tree <- tree
  }
  debug_log("asr_tree.f sigma2_used = %.6f", sigma2_used)
  
  # Node indices of the rescaled tree.
  Ntip_rescaled <- length(rescaled_tree$tip.label)
  Nnode_rescaled <- rescaled_tree$Nnode
  node_indices <- (Ntip_rescaled+1):(Ntip_rescaled+Nnode_rescaled)

  # ACE (REML, then ML) on the rescaled tree; when both fail, ACE on the original tree
  # (the model becomes BM), and finally fastAnc without CIs.
  ace_try <- function(tree_obj) {
    tryCatch(
      ace_fun(x = trait_vec, phy = tree_obj, type = "continuous", method = "REML"),
      error = function(e) tryCatch(
        ace_fun(x = trait_vec, phy = tree_obj, type = "continuous", method = "ML"),
        error = function(e2) NULL
      )
    )
  }
  ace_res <- ace_try(rescaled_tree)
  model_used <- chosen_model
  if (is.null(ace_res)) {
    message("ACE failed on rescaled tree; retrying on original tree.")
    ace_res <- ace_try(tree)
    if (is.null(ace_res)) {
      message("ACE failed on original tree; falling back to fastAnc without CIs.")
      anc_est <- tryCatch(
        phytools::fastAnc(rescaled_tree, trait_vec),
        error = function(e) NULL
      )
      if (is.null(anc_est)) {
        debug_log("asr_tree.f fastAnc failed")
      } else {
        debug_log("asr_tree.f fastAnc ok, n = %d", length(anc_est))
      }
      if (is.null(anc_est)) {
        stop("ACE and fastAnc failed.")
      }
      ace_res <- list(
        ace = anc_est,
        CI95 = matrix(NA_real_, nrow = length(anc_est), ncol = 2)
      )
    } else {
      rescaled_tree <- tree
      Ntip_rescaled <- length(rescaled_tree$tip.label)
      Nnode_rescaled <- rescaled_tree$Nnode
      node_indices <- (Ntip_rescaled+1):(Ntip_rescaled+Nnode_rescaled)
      model_used <- "BM"
    }
  }
  
  node_est <- ace_res$ace
  tip_est <- trait_vec[rescaled_tree$tip.label]
  debug_log("asr_tree.f node_est = %d, tip_est = %d", length(node_est), length(tip_est))
  
  # Node and tip tables. Node CIs are the analytical 95% CIs of ace().
  node_ci <- ace_res$CI95
  node_df <- data.frame(
    node = node_indices,
    estimate = as.numeric(node_est),
    ci_lower = as.numeric(node_ci[,1]),
    ci_upper = as.numeric(node_ci[,2]),
    stringsAsFactors = FALSE
  )
  # Tips carry no CI in the ace output: the observed value is used as both bounds.
  tip_df <- data.frame(
    tip = rescaled_tree$tip.label,
    estimate = as.numeric(tip_est[rescaled_tree$tip.label]),
    ci_lower = as.numeric(tip_est[rescaled_tree$tip.label]),
    ci_upper = as.numeric(tip_est[rescaled_tree$tip.label]),
    immediate_parent = NA,
    cumulative_derived = NA,
    direction = NA,
    parent_node_id = NA,
    parent_ci_lower = NA,
    parent_ci_upper = NA,
    stringsAsFactors = FALSE
  )
  
  # Derivedness tests: against the immediate parent and against all ancestors.
  # parent_map[i] is the parent node of node i (NA for the root).
  parent_map <- integer(Ntip_rescaled + Nnode_rescaled)
  parent_map[] <- NA_integer_
  for (j in seq_len(nrow(rescaled_tree$edge))) {
    parent_map[rescaled_tree$edge[j,2]] <- rescaled_tree$edge[j,1]
  }
  
  # A tip is derived when its value lies outside the CI of its parent node
  # (direction "up" or "down").
  for (i in seq_len(nrow(tip_df))) {
    tip_label <- tip_df$tip[i]
    tip_est_val <- tip_df$estimate[i]
    tip_idx <- which(rescaled_tree$tip.label == tip_label)
    parent_id <- parent_map[tip_idx]
    
    if (!is.na(parent_id)) {
      tip_df$parent_node_id[i] <- parent_id
      
      if (parent_id <= Ntip_rescaled) {
        parent_row <- tip_df[tip_df$tip == rescaled_tree$tip.label[parent_id], ]
      } else {
        parent_row <- node_df[node_df$node == parent_id, ]
      }
      parent_ci_low <- parent_row$ci_lower
      parent_ci_high <- parent_row$ci_upper
      
      if (!all(is.finite(c(parent_ci_low, parent_ci_high)))) next
      tip_df$immediate_parent[i] <- (tip_est_val > parent_ci_high) || (tip_est_val < parent_ci_low)
      
      if (tip_est_val > parent_ci_high) {
        tip_df$direction[i] <- "up"
      } else if (tip_est_val < parent_ci_low) {
        tip_df$direction[i] <- "down"
      } else {
        tip_df$direction[i] <- "none"
      }
      
      if (tip_df$immediate_parent[i]) {
        tip_df$parent_ci_lower[i] <- parent_ci_low
        tip_df$parent_ci_upper[i] <- parent_ci_high
      }
    }
    
    # Cumulatively derived: outside the CI of every ancestor node.
    anc_list <- integer(0)
    cur <- tip_idx
    while (!is.na(parent_map[cur])) {
      cur <- parent_map[cur]
      anc_list <- c(anc_list, cur)
    }
    anc_nodes <- anc_list[anc_list > Ntip_rescaled]
    
    if (length(anc_nodes) == 0) {
      tip_df$cumulative_derived[i] <- FALSE
    } else {
      anc_rows <- node_df[node_df$node %in% anc_nodes, ]
      if (nrow(anc_rows) == 0 || any(!is.finite(anc_rows$ci_lower)) || any(!is.finite(anc_rows$ci_upper))) {
        tip_df$cumulative_derived[i] <- NA
      } else {
        tip_df$cumulative_derived[i] <- all((tip_est_val > anc_rows$ci_upper) | (tip_est_val < anc_rows$ci_lower))
      }
    }
  }
  
  # Plot (ggtree). Branches are colored by the node estimate; tips carry the species
  # name, with the estimate in brackets when the tip is derived from its parent. Text and
  # symbol sizes shrink slowly with the number of tips, while the canvas size comes from
  # species_plot_profile()$asr.
  all_vals <- rbind(
    data.frame(node = seq_len(Ntip_rescaled), estimate = as.numeric(tip_est[rescaled_tree$tip.label])),
    data.frame(node = node_df$node, estimate = node_df$estimate)
  )
  tip_meta <- tip_df
  tip_meta$node <- seq_len(Ntip_rescaled)
  tip_meta$is_up <- !is.na(tip_meta$direction) & tip_meta$direction == "up"
  tip_meta$is_down <- !is.na(tip_meta$direction) & tip_meta$direction == "down"
  tip_meta$species_annotated <- ifelse(
    tip_meta$immediate_parent %in% TRUE,
    sprintf("plain(%s)~bold(\"[%.2f]\")", vapply(tip_meta$tip, deparse, character(1)), tip_meta$estimate),
    sprintf("plain(%s)", vapply(tip_meta$tip, deparse, character(1)))
  )

  # Parent nodes of derived tips get a label with the estimate and its CI. The labels are
  # joined into plot_df because geom_nodelab() needs the tree-attached data (x, y and isTip
  # from the %<+% join); a standalone data frame without those columns fails in its
  # internal `!isTip` check.
  parent_nodes_with_derived_tips <- unique(tip_df$parent_node_id[tip_df$immediate_parent %in% TRUE & !is.na(tip_df$parent_node_id)])
  node_df$node_label <- ifelse(
    node_df$node %in% parent_nodes_with_derived_tips,
    sprintf("%.2f (%.2f-%.2f)", node_df$estimate, node_df$ci_lower, node_df$ci_upper),
    NA_character_
  )

  plot_df <- dplyr::left_join(
    all_vals,
    tip_meta[, c("node", "is_up", "is_down", "cumulative_derived", "species_annotated")],
    by = "node"
  )
  plot_df$is_up[is.na(plot_df$is_up)] <- FALSE
  plot_df$is_down[is.na(plot_df$is_down)] <- FALSE
  plot_df$cumulative_derived[is.na(plot_df$cumulative_derived)] <- FALSE
  plot_df <- dplyr::left_join(plot_df, node_df[, c("node", "node_label")], by = "node")

  root_est_str <- if (all(is.finite(c(node_df$ci_lower[1], node_df$ci_upper[1])))) {
    sprintf("Root: %.2f (%.2f-%.2f)", node_df$estimate[1], node_df$ci_lower[1], node_df$ci_upper[1])
  } else {
    sprintf("Root: %.2f (CI unavailable; derivedness not assessed)", node_df$estimate[1])
  }

  # Text and symbol sizes decrease linearly with the number of tips, down to a floor,
  # so the labels do not overlap vertically in large trees.
  tip_lab_size   <- max(2.2, 4.5 - 0.008 * Ntip_rescaled)
  symbol_size    <- max(2.5, 5.0 - 0.008 * Ntip_rescaled)
  node_lab_size  <- max(2.0, 4.0 - 0.008 * Ntip_rescaled)

  vermillion_col <- as.character(tableau_10["vermillion"])
  # "green" rather than "bluish_green": the teal of the latter is hard to see against the
  # white background in the down-triangle.
  green_col <- as.character(tableau_10["green"])
  blue_col <- as.character(tableau_10["blue"])

  gg_base <- ggtree::ggtree(rescaled_tree, layout = "rectangular", color = "black", linewidth = 1.2) %<+% plot_df
  tree_depth <- max(gg_base$data$x, na.rm = TRUE)
  # The tip rows of the joined data position the symbol layers at explicit (x, y), because
  # the tiplab offset and position_nudge only move along the branch axis.
  tip_data <- gg_base$data[gg_base$data$isTip, ]

  p <- gg_base +
    ggtree::geom_tree(ggplot2::aes(color = estimate), linewidth = 0.7) +
    ggplot2::scale_color_gradient(
      low = "white", high = vermillion_col,
      name = paste0(trait, " (", model_used, ")"),
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(20, "lines"),
        barheight = ggplot2::unit(1.2, "lines"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ggtree::geom_tiplab(ggplot2::aes(label = species_annotated), parse = TRUE, size = tip_lab_size, offset = 0.08 * tree_depth, color = "black") +
    ggtree::theme_tree2() +
    # hexpand(direction = 1) pads the x-range on the label side only, which frees room on
    # the right for tip labels, values and symbols that would be clipped at the panel edge.
    ggtree::hexpand(0.28, direction = 1) +
    ggplot2::theme(
      legend.position = "none",  # the color bar is added below the caption with cowplot
      plot.margin = ggplot2::margin(10, 60, 10, 10),
      plot.subtitle = ggplot2::element_text(size = 16, color = "black", face = "bold"),
      plot.caption = ggplot2::element_text(size = 13, color = "black", hjust = 0.35),
      axis.text = ggplot2::element_text(size = 12, color = "black")
    ) +
    ggplot2::labs(
      subtitle = root_est_str,
      caption = paste0("Symbols: ▲ Tip Up / ▼ Tip Down (derived vs. immediate parent node), ",
                        "■ Cumulatively derived (derived vs. ALL ancestor nodes); ",
                        "parent node values shown for derived tips.")
    )

  # Derivation symbols sit on the tip row, to the right of the tip: the up/down triangle at
  # the first offset and the cumulative-derived square at the second offset when a triangle
  # is also present, at the first offset otherwise.
  sym_offset1 <- 0.025 * tree_depth
  sym_offset2 <- 0.052 * tree_depth
  has_tri <- tip_data$is_up | tip_data$is_down
  tip_data$sq_x <- ifelse(has_tri, tip_data$x + sym_offset2, tip_data$x + sym_offset1)

  if (any(tip_data$is_up)) {
    p <- p + ggplot2::geom_point(data = subset(tip_data, is_up),
                                  ggplot2::aes(x = x + sym_offset1, y = y),
                                  shape = 24, fill = vermillion_col, color = "black", size = symbol_size)
  }
  if (any(tip_data$is_down)) {
    p <- p + ggplot2::geom_point(data = subset(tip_data, is_down),
                                  ggplot2::aes(x = x + sym_offset1, y = y),
                                  shape = 25, fill = green_col, color = "black", size = symbol_size)
  }
  if (any(tip_data$cumulative_derived)) {
    p <- p + ggplot2::geom_point(data = subset(tip_data, cumulative_derived),
                                  ggplot2::aes(x = sq_x, y = y),
                                  shape = 22, fill = blue_col, color = "black", size = symbol_size)
  }
  if (any(!is.na(plot_df$node_label))) {
    p <- p + ggtree::geom_nodelab(ggplot2::aes(subset = !is.na(node_label), label = node_label),
                                   size = node_lab_size, fontface = "bold", hjust = -0.1, color = "black")
  }

  # The color bar is extracted from the plot and stacked below the whole panel, caption
  # included; ggplot2's own bottom legend would sit between the panel and the caption.
  legend_grob <- cowplot::get_legend(
    p + ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 14, color = "black", margin = ggplot2::margin(b = 6)),
      legend.text = ggplot2::element_text(size = 12, color = "black")
    )
  )
  p <- cowplot::plot_grid(p, legend_grob, ncol = 1, rel_heights = c(1, 0.06))

  list(
    chosen_model=chosen_model,
    node_table=node_df,
    tip_table=tip_df,
    sigma2_used=sigma2_used,
    ace=ace_res,
    rescaled_tree=rescaled_tree,
    plot_dims=asr_dims[c("width", "height")],
    plot=p
  )
}



contrast_plot_highlighted.f <- function(df, ordered_taxa) {
  plot_profile <- species_plot_profile(
    n_species = dplyr::n_distinct(df$species),
    n_taxa = dplyr::n_distinct(df$taxa)
  )
  contrast_dims <- plot_profile$contrast
  # Make distribution points larger and more legible
  base_pt_size <- max(3.5, contrast_dims$point_size * 1.5)
  h1_pt_size   <- base_pt_size * 1.6
  fop_pt_size  <- base_pt_size * 1.5
  global_median <- df$g_median[1]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = taxa)) +
    ggplot2::scale_y_discrete(limits = rev(ordered_taxa)) +
    ggplot2::geom_point(
      ggplot2::aes(color = taxa), shape = 16, size = base_pt_size, alpha = 0.35
    )

  if (!is.null(taxa_palette_values) && length(taxa_palette_values) > 0) {
    p <- p + ggplot2::scale_color_manual(values = taxa_palette_values, guide = "none")
  } else {
    p <- p + ggplot2::scale_color_discrete(guide = "none")
  }

  p <- p +
    ggnewscale::new_scale_color() +
    # Secondary FOP species: open circles (unfilled, colored border)
    ggplot2::geom_point(
      data = subset(df, is_secondary),
      ggplot2::aes(color = contrast_label),
      shape = 1, stroke = 1.3, size = fop_pt_size
    ) +
    # Canonical H1 species: filled circles (shape 21, solid fill + black stroke)
    ggplot2::geom_point(
      data = subset(df, is_canonical),
      ggplot2::aes(fill = contrast_label),
      shape = 21, color = "black", stroke = 1.2, size = h1_pt_size
    ) +
    ggplot2::scale_color_manual(
      name = "Secondary (FOP):",
      values = c("Top/High (FOP secondary)" = "red", "Bottom/Low (FOP secondary)" = "forestgreen"),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = "Canonical (H1):",
      values = c("Top/High (H1)" = "red", "Bottom/Low (H1)" = "forestgreen"),
      drop = FALSE
    ) +
    ggplot2::geom_vline(xintercept = global_median, linetype = "longdash", linewidth = 1.2, color = "grey30") +
    ggplot2::labs(
      x = capitalized_trait,
      y = paste(capitalized_taxon, "-", capitalized_clade),
      title = paste("Trait distribution with selected contrast pairs highlighted:", capitalized_trait),
      caption = "Global median: dashed line. Filled points: Canonical H1 pairs. Open circles: FOP secondary alternatives."
    ) +
    report_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = contrast_dims$title_size, hjust = 0, margin = ggplot2::margin(t = 20, b = 6), face = "bold"),
      axis.text.y = ggplot2::element_text(size = contrast_dims$axis_text_y),
      axis.title.x = ggplot2::element_text(size = contrast_dims$axis_title, margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(size = contrast_dims$axis_title, angle = 90, margin = ggplot2::margin(l = 20, r = 20)),
      plot.caption = ggplot2::element_text(size = contrast_dims$caption_size, hjust = 0, margin = ggplot2::margin(t = 20)),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  p
}


violin_highlighted.f <- function(df, ordered_taxa) {
  plot_profile <- species_plot_profile(
    n_species = dplyr::n_distinct(df$species),
    n_taxa = dplyr::n_distinct(df$taxa)
  )
  violin_dims <- plot_profile$violin
  base_pt_size <- max(2.2, violin_dims$point_size * 1.3)
  h1_pt_size   <- base_pt_size * 1.6
  fop_pt_size  <- base_pt_size * 1.5
  median_value <- median(df$value, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = taxa)) +
    ggplot2::scale_y_discrete(limits = rev(ordered_taxa))

  # Violin fill with family palette
  if (!is.null(dark_palette_values) && length(dark_palette_values) > 0) {
    p <- p +
      ggplot2::geom_violin(
        ggplot2::aes(fill = taxa), color = "black", width = 0.5, trim = TRUE, scale = "width", alpha = 0.20
      ) +
      ggplot2::scale_fill_manual(values = dark_palette_values, guide = "none") +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(
        ggplot2::aes(color = taxa), shape = 16, size = base_pt_size, alpha = 0.35,
        position = ggplot2::position_jitter(height = violin_dims$jitter_height)
      ) +
      ggplot2::scale_color_manual(values = taxa_palette_values, guide = "none")
  } else {
    p <- p +
      ggplot2::geom_violin(
        ggplot2::aes(fill = taxa), color = "black", width = 0.5, trim = TRUE, scale = "width", alpha = 0.20
      ) +
      ggplot2::scale_fill_discrete(guide = "none") +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(
        ggplot2::aes(color = taxa), shape = 16, size = base_pt_size, alpha = 0.35,
        position = ggplot2::position_jitter(height = violin_dims$jitter_height)
      ) +
      ggplot2::scale_color_discrete(guide = "none")
  }

  p <- p +
    ggnewscale::new_scale_color() +
    # Secondary FOP species: open circles (shape 1, unfilled, jittered)
    ggplot2::geom_point(
      data = subset(df, is_secondary),
      ggplot2::aes(color = contrast_label),
      shape = 1, stroke = 1.2, size = fop_pt_size,
      position = ggplot2::position_jitter(height = violin_dims$jitter_height)
    ) +
    # Canonical H1 species: filled circles (shape 21, jittered)
    ggplot2::geom_point(
      data = subset(df, is_canonical),
      ggplot2::aes(fill = contrast_label),
      shape = 21, color = "black", stroke = 1.1, size = h1_pt_size,
      position = ggplot2::position_jitter(height = violin_dims$jitter_height)
    ) +
    ggplot2::scale_color_manual(
      name = "Secondary (FOP):",
      values = c("Top/High (FOP secondary)" = "red", "Bottom/Low (FOP secondary)" = "forestgreen"),
      drop = FALSE,
      guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 2)
    ) +
    ggplot2::scale_fill_manual(
      name = "Canonical (H1):",
      values = c("Top/High (H1)" = "red", "Bottom/Low (H1)" = "forestgreen"),
      drop = FALSE,
      guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1)
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = median_value), linewidth = 1, color = "grey30", alpha = 0.5) +
    ggplot2::labs(x = capitalized_trait) +
    report_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = violin_dims$axis_text_y, hjust = 0.5, face = "bold"),
      axis.title.x = ggplot2::element_text(size = violin_dims$axis_title, face = "bold", margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.15, "cm")
    )
  p
}
