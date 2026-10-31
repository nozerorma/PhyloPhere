# Shared visual defaults for the curated reports. HTML uses the existing
# Helvetica Neue / Helvetica / Arial / sans-serif stack; plots use its portable
# sans-serif equivalent. Individual figures may adjust sizes for their content.
report_palette <- function(n, palette = "ggthemes::Tableau_10") {
  base <- as.character(paletteer::paletteer_d(palette))
  if (!n) return(character())
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

report_theme <- function(base_size = 13, base_family = "sans", ...) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family, ...) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 20, r = 20)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(colour = "grey40", hjust = 0, margin = ggplot2::margin(b = 10)),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

set_report_style <- function() {
  # Cover discrete plots without an explicit semantic/taxon colour mapping.
  options(
    ggplot2.discrete.colour = function(...) ggplot2::discrete_scale("colour", palette = report_palette, ...),
    ggplot2.discrete.fill = function(...) ggplot2::discrete_scale("fill", palette = report_palette, ...)
  )
  ggplot2::theme_set(report_theme())
  invisible(NULL)
}
