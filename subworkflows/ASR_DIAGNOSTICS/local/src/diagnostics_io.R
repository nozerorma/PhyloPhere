# Input contract for report 9. Read only report columns; never load domain tallies.
asr_fail <- function(message) stop(paste0('[asr_diagnostics] ', message), call. = FALSE)
asr_types <- c('convergent', 'divergent', 'single', 'no_change',
               'symplesiomorphic', 'cross_hypothesis')

asr_number <- function(x, name, optional = FALSE, integer = FALSE, lower = 0, upper = Inf) {
  value <- suppressWarnings(as.numeric(x))
  missing <- is.na(x) | x == ''
  bad <- (!missing & (!is.finite(value) | value < lower | value > upper)) |
    (missing & !optional)
  if (integer) bad <- bad | (!missing & is.finite(value) & value != floor(value))
  if (any(bad)) asr_fail(paste('Invalid', name))
  value[missing] <- NA_real_
  value
}

asr_coverage <- function(df) {
  positions <- split(df$msa_pos, df$gene)
  lapply(positions, function(x) {
    x <- sort(unique(x))
    text <- paste0(paste(sprintf('%.0f', x), collapse = '\n'), '\n')
    list(positions = length(x), sha256 = digest::digest(text, algo = 'sha256', serialize = FALSE))
  })
}

asr_check_coverage <- function(actual, expected) {
  if (!setequal(names(actual), names(expected))) asr_fail('Master gene coverage differs from completion record')
  for (gene in names(actual)) {
    a <- actual[[gene]]; b <- expected[[gene]]
    if (!identical(a$sha256, b$sha256) || !identical(as.numeric(a$positions), as.numeric(b$positions)))
      asr_fail(paste('Master position coverage differs for', gene))
  }
}

load_asr_diagnostics <- function(directory) {
  master <- file.path(directory, 'caas_convergence_master.csv')
  stamp <- file.path(directory, 'observed_scoring.json')
  if (!file.exists(master) || !file.exists(stamp))
    asr_fail('Requires caas_convergence_master.csv and observed_scoring.json; regenerate legacy observed results')
  record <- jsonlite::read_json(stamp, simplifyVector = FALSE)
  if (!identical(record$version, 1L) || !identical(record$complete, TRUE) || is.null(record$genes))
    asr_fail('Invalid observed completion record')
  batches <- record$batches
  if (is.null(batches)) asr_fail('Missing observed batch provenance')
  settings <- c('model', 'posterior_threshold', 'tree_sha256', 'taxid_sha256', 'weights_sha256', 'weighted')
  genes <- list()
  excluded <- unlist(record$excluded_genes, use.names = FALSE)
  for (batch in batches) {
    if (!identical(batch$version, 1L) || !identical(batch$complete, TRUE) ||
        !identical(batch$score_rule, 'caas_pos_species_pss_v4') || is.null(batch$genes))
      asr_fail('Incomplete or incompatible observed batch')
    if (is.null(batch$model) || is.null(batch$tree_sha256) || is.null(batch$weighted) ||
        length(batch$posterior_threshold) != 1 || !is.numeric(batch$posterior_threshold))
      asr_fail('Missing scientific settings in observed batch')
    asr_number(batch$posterior_threshold, 'recorded posterior threshold', upper = 1)
    if (any(names(batch$genes) %in% names(genes))) asr_fail('Overlapping genes in observed batches')
    if (any(!vapply(settings, function(key) identical(batch[[key]], batches[[1]][[key]]), logical(1))))
      asr_fail('Inconsistent scientific settings between observed batches')
    genes <- c(genes, batch$genes)
    excluded <- c(excluded, unlist(batch$excluded_genes, use.names = FALSE))
  }
  if (any(names(genes) %in% excluded)) asr_fail('Genes are both scored and excluded')
  asr_check_coverage(genes, record$genes)
  if (!length(batches) && length(record$genes)) asr_fail('Nonempty results without batch provenance')

  header <- names(readr::read_csv(master, n_max = 0, name_repair = 'minimal', show_col_types = FALSE))
  required <- c('gene', 'msa_pos', 'side', 'asr_path_score')
  optional <- c('derived_agreement', 'convergence_type', 'agreement_ambiguous', 'n_hypotheses')
  if (anyDuplicated(header) || !all(required %in% header)) asr_fail('Invalid master header')
  df <- readr::read_csv(master, col_select = tidyselect::all_of(intersect(header, c(required, optional))),
                        col_types = readr::cols(.default = readr::col_character()),
                        name_repair = 'minimal', num_threads = 1, progress = FALSE)
  if (nrow(readr::problems(df))) asr_fail('Malformed master CSV')
  for (key in setdiff(optional, names(df))) df[[key]] <- rep(NA_character_, nrow(df))
  if (any(is.na(df$gene) | trimws(df$gene) == '') || any(!df$side %in% c('top', 'bottom', 'none')))
    asr_fail('Missing gene or invalid side')
  df$msa_pos <- asr_number(df$msa_pos, 'msa_pos', integer = TRUE)
  df$asr_path_score <- asr_number(df$asr_path_score, 'asr_path_score', upper = 1)
  df$derived_agreement <- asr_number(df$derived_agreement, 'derived_agreement', optional = TRUE, upper = 1)
  df$n_hypotheses <- asr_number(df$n_hypotheses, 'n_hypotheses', optional = TRUE, integer = TRUE, lower = 1)
  if (any(!is.na(df$convergence_type) & !df$convergence_type %in% asr_types)) asr_fail('Invalid convergence_type')
  flag <- tolower(df$agreement_ambiguous)
  if (any(!is.na(flag) & !flag %in% c('true', 'false'))) asr_fail('Invalid agreement_ambiguous')
  df$agreement_ambiguous <- ifelse(is.na(flag), NA, flag == 'true')
  if (anyDuplicated(df[c('gene', 'msa_pos', 'side')])) asr_fail('Duplicate gene/position/side')
  conflicts <- df %>% dplyr::group_by(gene, msa_pos) %>%
    dplyr::summarise(conflict = any(side == 'none') && dplyr::n() > 1, .groups = 'drop')
  if (any(conflicts$conflict) || any(df$side == 'none' & df$asr_path_score != 0))
    asr_fail('Conflicting none-side rows')
  asr_check_coverage(asr_coverage(df), record$genes)
  list(data = df, record = record,
       posterior_threshold = if (length(batches)) batches[[1]]$posterior_threshold else NA_real_,
       model = if (length(batches)) batches[[1]]$model else 'not recorded (empty result)',
       excluded_genes = length(unique(excluded)))
}

asr_percent <- function(x) if (all(is.na(x))) NA_real_ else 100 * mean(x, na.rm = TRUE)
