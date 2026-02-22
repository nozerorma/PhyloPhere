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
