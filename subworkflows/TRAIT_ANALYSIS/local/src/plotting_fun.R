# plotting_fun.R — Palette lookup, tree-annotation helpers and figure-size profiles.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
# =============================================================================
# Sourced by: commons.R (itself sourced by the trait-analysis Rmd reports)
#
# Defines:
#   get_palette_values(), resolve_taxa_palette()  palette lookup by name
#   clamp_value(), compute_ring_axis_limits()     small numeric helpers
#   branch_trait_node_values()                    per-node values of a branch trait
#   fan_label_offsets()                           label placement on fan trees
#   species_plot_profile()                        figure sizes as a function of
#                                                 the number of species and taxa
# =============================================================================

# ── Palettes ──────────────────────────────────────────────────────────────────

# Return the palette object called `palette_name` (default: the global
# `color_palette`), or NULL when it is not defined. The lookup happens at call
# time, so palettes.R may be sourced before or after this file.
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

# Named color vector (one color per distinct taxon in `taxa_values`). Taxa found
# in the clade palette keep its color; the rest get a color from the fallback
# palette chosen by the sum of the character codes of the taxon name, so the same
# taxon receives the same color in every plot.
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

# ── Numeric helpers ───────────────────────────────────────────────────────────

# Restrict `x` to the interval [lower, upper].
clamp_value <- function(x, lower, upper) {
  max(lower, min(upper, x))
}

# Axis limits c(0, upper) of a ring plot: the largest finite value plus `headroom`
# (a multiplicative margin), never below `min_limit`.
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

# ── Branch trait on a tree ────────────────────────────────────────────────────

# Per-node values of a branch-colour trait, for the branch colors of the fan plot.
#
# Trait tables routinely break a direct fastAnc() call on the plotting data:
#   * duplicated species rows give a vector longer than the number of tips, and
#     ace() stops with "length of phenotypic and of phylogenetic data do not match";
#   * a branch trait is usually recorded for a subset of the species carrying the
#     primary trait, while ancestral reconstruction needs complete data over the
#     tree it is given.
#
# The function collapses the data to one finite value per species, reconstructs
# on the subtree that has data, and maps each subtree node onto the full plotting
# tree through the clade it subtends. Nodes of the full tree that subtend fewer
# than two covered species get no value (the plot draws them in the `na.value`
# color), because there is nothing to reconstruct there.
#
# Returns a data.frame(node, BR) over the tips and internal nodes of `tree`, or
# NULL when the trait cannot support a reconstruction at all.
branch_trait_node_values <- function(tree, species, values) {
  df <- data.frame(
    species = as.character(species),
    value = suppressWarnings(as.numeric(values)),
    stringsAsFactors = FALSE
  )
  n_raw <- nrow(df)
  df <- df[df$species %in% tree$tip.label & is.finite(df$value), , drop = FALSE]
  if (nrow(df) == 0) {
    debug_log("branch_trait_node_values: no finite values on tree tips (from %d rows)", n_raw)
    return(NULL)
  }

  # One value per species; duplicated rows are averaged rather than dropped
  # arbitrarily, so a conflicting duplicate does not silently pick a side.
  df <- stats::aggregate(value ~ species, data = df, FUN = mean)
  if (nrow(df) < 3) {
    debug_log("branch_trait_node_values: only %d species with data; skipping ASR", nrow(df))
    return(NULL)
  }

  n_tip_full <- length(tree$tip.label)
  sub_tree <- if (nrow(df) < n_tip_full) {
    ape::drop.tip(tree, setdiff(tree$tip.label, df$species))
  } else {
    tree
  }
  trait_vec <- stats::setNames(df$value, df$species)[sub_tree$tip.label]
  debug_log("branch_trait_node_values: %d/%d tips covered (%d raw rows)",
            length(trait_vec), n_tip_full, n_raw)

  anc <- tryCatch(
    phytools::fastAnc(sub_tree, trait_vec),
    error = function(e) {
      debug_log("branch_trait_node_values: fastAnc failed (%s)", conditionMessage(e))
      NULL
    }
  )
  if (is.null(anc)) return(NULL)

  tip_rows <- data.frame(
    node = match(names(trait_vec), tree$tip.label),
    BR = unname(trait_vec)
  )

  n_tip_sub <- length(sub_tree$tip.label)
  sub_nodes <- seq.int(n_tip_sub + 1L, n_tip_sub + sub_tree$Nnode)
  node_rows <- lapply(sub_nodes, function(nd) {
    desc <- phytools::getDescendants(sub_tree, nd)
    clade_tips <- sub_tree$tip.label[desc[desc <= n_tip_sub]]
    full_node <- ape::getMRCA(tree, clade_tips)
    est <- anc[[as.character(nd)]]
    if (is.null(full_node) || is.null(est)) return(NULL)
    data.frame(node = as.integer(full_node), BR = as.numeric(est))
  })
  node_rows <- do.call(rbind, node_rows[!vapply(node_rows, is.null, logical(1))])

  out <- rbind(tip_rows, node_rows)
  out <- out[!duplicated(out$node), , drop = FALSE]
  out$node <- as.numeric(out$node)
  debug_log("branch_trait_node_values: %d of %d tree nodes assigned a value",
            nrow(out), n_tip_full + tree$Nnode)
  out
}

# ── Fan-plot label placement ──────────────────────────────────────────────────

# Radial offsets for the taxon labels and phylopic images of a fan plot, placed
# just outside the stack of rings.
#
# The offsets are measured on the built plot rather than tabulated, because the
# two quantities to align use different units:
#   * geom_fruit sizes each ring as a fraction of the tree's plotted x-range and
#     stacks the rings outward, so the outer edge depends on every ring added;
#   * geom_cladelab takes an absolute offset from the tree's own x-range, which
#     does not know about the rings outside it.
# The plotted x-range is whatever ggtree produces and need not match the
# branch-length depth, so a fixed offset suits one dataset only.
#
# The function builds the ring plot, reads where the rings end, and returns
# offsets that put the labels (`text`) and then the images (`phylopic`) outside.
#
# `ring_plot`      ggplot of the fan tree with its rings
# `labels`         taxon labels (the longest sets the radial room for the text)
# `label_fontsize` label font size in ggplot units (mm)
# `canvas_in`      side of the panel in inches (the fan is bounded by plot
#                  height), used to convert the label length into x-units
# `ring_gap`, `image_gap`  gaps before the text and before the image, as
#                  fractions of the tree radius
# `image_size`     phylopic diameter as a fraction of the canvas
# Returns list(tip_radius, outer_radius, text, phylopic).
fan_label_offsets <- function(ring_plot, labels = character(), label_fontsize = 8,
                              canvas_in = 20, ring_gap = 0.10, image_gap = 0.06,
                              image_size = 0.04) {
  tip_radius <- suppressWarnings(max(ring_plot$data$x, na.rm = TRUE))
  if (!is.finite(tip_radius) || tip_radius <= 0) tip_radius <- 1

  built <- ggplot2::ggplot_build(ring_plot)
  ring_x <- unlist(lapply(built$data, function(d) c(d$x, d$xmax)), use.names = FALSE)
  ring_x <- suppressWarnings(as.numeric(ring_x))
  ring_x <- ring_x[is.finite(ring_x)]
  outer_radius <- if (length(ring_x)) max(ring_x) else tip_radius

  # Distance from the cladelab origin (the tree tips) to the outer edge of the last ring.
  ring_span <- max(0, outer_radius - tip_radius)

  text_offset <- ring_span + ring_gap * tip_radius

  # Radial room the taxon names need before the phylopics start. ggplot font
  # sizes are in mm and a character advances about 0.55 of the font height; the
  # panel spans `outer_radius` x-units over half the canvas, which converts
  # inches to x-units.
  units_per_inch <- outer_radius / max(canvas_in / 2, 1e-6)
  label_chars <- if (length(labels)) {
    max(nchar(as.character(labels)), na.rm = TRUE)
  } else {
    8
  }
  label_span <- (label_chars * 0.55 * label_fontsize / 25.4) * units_per_inch

  # geom_cladelab grows its labels outward from the offset, so the phylopic must
  # clear the whole longest name plus half its own diameter (the image is
  # centered on its offset).
  image_radius <- (image_size * canvas_in / 2) * units_per_inch
  image_offset <- text_offset + label_span + image_radius + image_gap * tip_radius

  debug_log(paste("fan_label_offsets: tip_radius=%.2f outer_radius=%.2f",
                  "text_offset=%.2f image_offset=%.2f (longest label %d chars)"),
            tip_radius, outer_radius, text_offset, image_offset, label_chars)

  list(
    tip_radius = tip_radius,
    outer_radius = outer_radius,
    text = text_offset,
    phylopic = image_offset
  )
}

# ── Figure-size profiles ──────────────────────────────────────────────────────

# Figure dimensions and text/point sizes as a function of dataset size, so plots
# stay legible from ~30 to several hundred species.
#
# `n_species`, `n_taxa` number of tips and of taxa (e.g. families) shown;
# `n_rings` number of annotation rings around the fan tree.
# Returns list(contrast, violin, asr, tree), one named list of settings per plot
# type (sizes in inches unless the name says otherwise), plus the three counts.
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

  # Fan canvas relative to a ~30-species reference (tree_height about 13.5 in),
  # and how tightly the taxon labels are packed around the annotation ring.
  tree_scale <- tree_height / 13.5
  taxa_crowding <- clamp_value(1.15 - 0.014 * n_taxa, 0.72, 1.0)

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
      # Width grows faster than height (unlike the other profiles, which are
      # taller than wide): species-name labels, branch structure, node-value text
      # (e.g. "0.63 (0.38-0.86)") and dense clusters of derived tips all compete
      # for horizontal space as the species count rises.
      width = clamp_value(14 + 0.14 * n_species, 14, 42),
      height = clamp_value(12 + 0.10 * n_species, 12, 32),
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

      # --- Ring geometry (fractions of the tree's plotted x-range) ---------
      # Radial thickness is not what crowds as species are added (the sectors
      # get thinner angularly while each ring keeps its radial room), so these
      # depend on the number of rings to stack, not on n_species.
      fruit_pwidth = clamp_value(0.62 - 0.05 * pmax(0, n_rings - 1), 0.35, 0.62),
      fruit_gap = 0.03,
      secondary_offset = 0.03,
      taxa_bar_offset = 0.15,
      taxa_bar_pwidth = 0.09,

      # --- Annotation sizing ------------------------------------------------
      # Text, legend keys and phylopics are drawn at absolute sizes (pt/mm/cm)
      # on a canvas whose side grows with n_species (tree_height above), so they
      # scale up with the canvas. What crowds the annotation ring is the number
      # of taxa (one label each), hence the mild taxon penalty.
      axis_text_size = clamp_value(6.5 * tree_scale, 5.0, 14.0),
      axis_nbreak = if (n_species > 80) 1 else 2,
      branch_label_size = clamp_value(7.6 * tree_scale * taxa_crowding, 5.0, 22.0),
      image_size = clamp_value(0.042 * taxa_crowding, 0.022, 0.050),
      legend_title_size = clamp_value(14.0 * tree_scale, 12.0, 30.0),
      legend_text_size = clamp_value(12.0 * tree_scale, 10.0, 26.0),
      legend_key_cm = clamp_value(1.0 * tree_scale, 1.0, 2.6),
      legend_spacing_cm = clamp_value(1.5 * tree_scale, 1.0, 3.0),
      caption_size = clamp_value(13.0 * tree_scale, 11.0, 26.0),

      # Radial gaps between the outermost ring and the taxon text, and between
      # that text and the phylopic, as fractions of the tree radius (the
      # ring_gap and image_gap of fan_label_offsets()).
      label_ring_gap = 0.10,
      label_image_gap = 0.06,

      # Per-species annotations tighten as species are added, because the arc
      # available to each tip shrinks faster than the canvas grows.
      n_label_size = clamp_value(6.0 - 0.020 * n_species, 2.3, 6.0),
      n_label_nudge = clamp_value(8.4 + 0.60 * pmax(0, n_rings - 2) - 0.040 * n_species, 3.2, 9.8),
      asterisk_size = clamp_value(7.0 - 0.025 * n_species, 2.5, 7.0),
      asterisk_nudge = clamp_value(7.2 + 0.55 * pmax(0, n_rings - 2) - 0.032 * n_species, 2.8, 8.6),
      asterisk_y = clamp_value(0.18 - 0.0010 * n_species, 0.07, 0.18)
    )
  )
}
