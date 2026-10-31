# Strict, stable report tables: identifiers and residue strings remain literal.
pp_table <- function(path, required = character()) {
  if (is.null(path) || !file.exists(path)) stop('Missing report input: ', path)
  df <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()),
                        na = '', name_repair = 'minimal', show_col_types = FALSE,
                        num_threads = 1, progress = FALSE)
  if (anyDuplicated(names(df)) || !all(required %in% names(df)) || nrow(readr::problems(df)))
    stop('Malformed report table: ', path)
  df
}
pp_numeric <- function(df, columns) {
  for (col in intersect(columns, names(df))) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  df
}
pp_discovery <- function(path) {
  df <- pp_table(path, c('Gene', 'Position', 'side', 'asr_path_score'))
  df <- pp_numeric(df, c('Position', 'asr_path_score'))
  if (!'convergence_type' %in% names(df)) df$convergence_type <- rep('unknown', nrow(df))
  df$convergence_type[is.na(df$convergence_type)] <- 'unknown'
  if (anyDuplicated(df[c('Gene','Position','side')]) ||
      any(!is.finite(df$Position) | df$Position < 0 | df$Position != floor(df$Position)) ||
      any(!is.finite(df$asr_path_score) | df$asr_path_score < 0 | df$asr_path_score > 1) ||
      any(!df$side %in% c('top','bottom','none'))) stop('Invalid discovery keys/scores')
  df
}
pp_plot <- function(plot, filename, width = 10, height = 6) {
  print(plot)
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = 180)
}
pp_show <- function(df, caption) {
  if (!nrow(df)) { cat('No records.\n\n'); return(invisible(NULL)) }
  print(knitr::kable(df, caption = caption) %>% kableExtra::kable_styling(
    bootstrap_options = c('striped','hover','condensed')))
}
