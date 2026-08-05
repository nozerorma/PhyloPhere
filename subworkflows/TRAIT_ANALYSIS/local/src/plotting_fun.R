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

# Branch-colour trait -> per-node values for the fan plot.
#
# The naive version of this (build a named vector straight off the plotting data
# frame and hand it to fastAnc) breaks on two things that real trait tables hit
# routinely:
#   * duplicated species rows -- the vector then has more entries than the tree
#     has tips and ace() dies with "length of phenotypic and of phylogenetic data
#     do not match";
#   * partial coverage -- a branch trait is usually recorded for a subset of the
#     species carrying the primary trait (LQ is known for ~60% of primates), and
#     ASR needs complete data over whatever tree it is given.
#
# So: collapse to one finite value per species, reconstruct on the subtree that
# actually has data, then carry each subtree node's estimate back onto the full
# plotting tree by matching the clades they share. Full-tree nodes subtending
# fewer than two covered species get no value and render in `na.value` grey,
# which is the honest depiction -- there is nothing to reconstruct there.
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
  # arbitrarily, so a genuinely conflicting duplicate does not silently pick a side.
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

# Where to put the family text and phylopic rings of a fan plot.
#
# This has to be measured rather than tabulated, because the two things being
# lined up are quoted in different units:
#   * geom_fruit sizes each ring as a FRACTION of the tree's plotted x-range and
#     stacks them outward, so the outer edge depends on every ring added;
#   * geom_cladelab takes an ABSOLUTE offset, applied from the tree's own
#     x-range -- which knows nothing about the rings sitting outside it.
# Worse, the tree's plotted x-range is not the branch-length depth: ggtree draws
# this fan from a treedata object without usable branch lengths, so the range is
# a node-depth count (~13 for 157 primates, ~9 for 31) that shifts with tree
# shape. Any hard-coded offset is therefore tuned to one dataset and wrong for
# the next -- which is exactly how the labels ended up marooned far outside the
# rings on the diet trees while sitting on top of them on the cancer trees.
#
# So: build the ring stack, read off where it really ends, and return offsets
# that drop the labels just outside it.
#
# `canvas_in` is the panel's side in inches (the fan is bounded by plot height),
# used to convert the label's physical length into plot x-units.
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

  # Distance from the cladelab origin (the tree tips) out past the last ring.
  ring_span <- max(0, outer_radius - tip_radius)

  text_offset <- ring_span + ring_gap * tip_radius

  # Radial room the family names themselves need before the phylopics start.
  # ggplot font sizes are in mm; a character advances roughly 0.55 of the font
  # height. The panel spans `outer_radius` x-units over half the canvas, which
  # converts inches to x-units.
  units_per_inch <- outer_radius / max(canvas_in / 2, 1e-6)
  label_chars <- if (length(labels)) {
    max(nchar(as.character(labels)), na.rm = TRUE)
  } else {
    8
  }
  label_span <- (label_chars * 0.55 * label_fontsize / 25.4) * units_per_inch

  # geom_cladelab grows its labels outward from the offset, so the phylopic has
  # to clear the whole length of the longest name, plus half its own diameter
  # (the image is centred on its offset, not anchored at its inner edge).
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

  # How much bigger the fan canvas is than the ~30-species reference case the
  # original absolute sizes were tuned against (tree_height 13.4in there), and
  # how tightly the family labels are packed around the annotation ring.
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
      # width's old 24-inch cap was reached by any tree above ~125 species
      # (14 + 0.08*125 = 24) and never grew past it -- for a full primate
      # dataset (150-200+ tips) that produced a canvas TALLER than it was
      # wide (height's cap of 28 sits comfortably above it), even though this
      # plot needs MORE horizontal room than height as species count grows:
      # long species-name labels, branch structure, AND node-value annotation
      # text (e.g. "0.63 (0.38-0.86)") all compete for the same horizontal
      # space, and dense derived-tip clusters need width to avoid overlapping
      # each other. Raised coefficient and cap so width keeps scaling well
      # past where it used to plateau, and grows faster than height (unlike
      # every other profile above, which are taller-than-wide by design).
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
      # These used to taper off with n_species, which is why a large tree ended
      # up with a thin phenotype ring pushed far away from an equally thin taxon
      # ring: at 157 species fruit_pwidth bottomed out at 0.24 while
      # taxa_bar_offset maxed out at 0.26, so the GAP between the two rings was
      # wider than the phenotype ring itself. Radial thickness is not what gets
      # crowded as species are added -- the sectors get thinner angularly, and
      # the ring keeps exactly the same radial room -- so these are now fixed by
      # how many rings have to be stacked, not by how many species there are.
      fruit_pwidth = clamp_value(0.62 - 0.05 * pmax(0, n_rings - 1), 0.35, 0.62),
      fruit_gap = 0.03,
      secondary_offset = 0.03,
      taxa_bar_offset = 0.15,
      taxa_bar_pwidth = 0.09,

      # --- Annotation sizing ------------------------------------------------
      # Text, legend keys and phylopics are drawn at absolute sizes (pt/mm/cm)
      # onto a canvas whose side grows with n_species (tree_height above). The
      # previous profile shrank them as n_species grew, so large trees were
      # penalised twice over -- a 4.5pt family label on the 30in canvas of a
      # 157-species tree renders about a third the apparent size of the 7.6pt
      # label on the 17in canvas of a 31-species tree, which is precisely why
      # the settings tuned on small trait sets fell apart on large ones.
      #
      # What actually crowds the annotation ring is the number of TAXA (one
      # label each), and the ring's circumference grows with the canvas anyway,
      # so these scale UP with the canvas and take only a mild taxon penalty.
      axis_text_size = clamp_value(6.5 * tree_scale, 5.0, 14.0),
      axis_nbreak = if (n_species > 80) 1 else 2,
      branch_label_size = clamp_value(7.6 * tree_scale * taxa_crowding, 5.0, 22.0),
      image_size = clamp_value(0.042 * taxa_crowding, 0.022, 0.050),
      legend_title_size = clamp_value(14.0 * tree_scale, 12.0, 30.0),
      legend_text_size = clamp_value(12.0 * tree_scale, 10.0, 26.0),
      legend_key_cm = clamp_value(1.0 * tree_scale, 1.0, 2.6),
      legend_spacing_cm = clamp_value(1.5 * tree_scale, 1.0, 3.0),
      caption_size = clamp_value(13.0 * tree_scale, 11.0, 26.0),

      # Radial breathing room between the outermost ring and the family text,
      # and between that text and the phylopic, as fractions of the tree radius.
      # Consumed by fan_label_offsets(), which turns them into the absolute
      # offsets geom_cladelab wants.
      label_ring_gap = 0.10,
      label_image_gap = 0.06,

      # Per-species annotations: unlike the family labels these DO get tighter
      # as species are added, because the arc available to each tip shrinks
      # faster than the canvas grows, so they keep an n_species penalty.
      n_label_size = clamp_value(6.0 - 0.020 * n_species, 2.3, 6.0),
      n_label_nudge = clamp_value(8.4 + 0.60 * pmax(0, n_rings - 2) - 0.040 * n_species, 3.2, 9.8),
      asterisk_size = clamp_value(7.0 - 0.025 * n_species, 2.5, 7.0),
      asterisk_nudge = clamp_value(7.2 + 0.55 * pmax(0, n_rings - 2) - 0.032 * n_species, 2.8, 8.6),
      asterisk_y = clamp_value(0.18 - 0.0010 * n_species, 0.07, 0.18)
    )
  )
}
