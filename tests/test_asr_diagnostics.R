# Validate the consumer contract and render real report HTML on synthetic results.
suppressPackageStartupMessages(library(dplyr))
root <- normalizePath('.')
source('subworkflows/ASR_DIAGNOSTICS/local/src/diagnostics_io.R')
scratch <- commandArgs(trailingOnly = TRUE)[1]
dir.create(scratch, recursive = TRUE)
base <- data.frame(gene = c('001.2', '001.2', 'G.3', 'G.4', 'G.5'), msa_pos = c(0, 0, 1, 2, 3),
  side = c('top', 'bottom', 'none', 'none', 'none'), asr_path_score = c(.8, .3, 0, 0, 0),
  derived_agreement = c(1, 1, NA, 1, 1), convergence_type = c('convergent', 'convergent',
  'no_change', 'symplesiomorphic', 'cross_hypothesis'), agreement_ambiguous = c(FALSE, TRUE, FALSE, FALSE, FALSE),
  n_hypotheses = 2)
record_for <- function(df) {
  genes <- asr_coverage(df)
  # jsonlite encodes the empty named mapping as an object, as Python producers do.
  if (!length(genes)) genes <- structure(list(), names = character())
  batch <- list(version = 1L, complete = TRUE, score_rule = 'caas_pos_species_pss_v4', genes = genes,
    excluded_genes = list('excluded'), model = 'lg', posterior_threshold = .2,
    tree_sha256 = 'tree', weighted = TRUE, weights_sha256 = 'weights', taxid_sha256 = NULL)
  list(version = 1L, complete = TRUE, genes = genes, excluded_genes = list('excluded'), batches = list(batch))
}
write_case <- function(name, df = base, record = record_for(df)) {
  d <- file.path(scratch, name); dir.create(d, recursive = TRUE)
  write.csv(df, file.path(d, 'caas_convergence_master.csv'), row.names = FALSE, na = '')
  jsonlite::write_json(record, file.path(d, 'observed_scoring.json'), auto_unbox = TRUE, null = 'null')
  d
}
fails <- function(name, df = base, record = record_for(df), message) {
  d <- write_case(name, df, record)
  error <- tryCatch({load_asr_diagnostics(d); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(message, conditionMessage(error), fixed = TRUE))
  cat('PASS rejected:', name, '\n')
}
good <- load_asr_diagnostics(write_case('valid'))
stopifnot(identical(good$data$gene, base$gene), good$posterior_threshold == .2,
          good$excluded_genes == 1, nrow(good$data) == 5)
stopifnot(identical(asr_coverage(base)$`001.2`$sha256,
  digest::digest('0\n', algo='sha256', serialize=FALSE)))
for (value in c('Inf', '-Inf', 'NaN', 'oops', '-0.1', '1.1', '')) {
  df <- base; df$asr_path_score[1] <- value
  fails(paste0('score_', match(value, c('Inf','-Inf','NaN','oops','-0.1','1.1',''))), df, message = 'Invalid asr_path_score')
}
df <- base; df$msa_pos[1] <- .5; fails('fractional_position', df, message='Invalid msa_pos')
df <- base; df$side[1] <- 'bad'; fails('side', df, message='invalid side')
df <- base; df$agreement_ambiguous[1] <- 'yes'; fails('flag', df, message='Invalid agreement_ambiguous')
df <- base; df$derived_agreement[1] <- 1.2; fails('agreement', df, message='Invalid derived_agreement')
df <- base; df$n_hypotheses[1] <- 0; fails('hypotheses', df, message='Invalid n_hypotheses')
df <- base; df$convergence_type[1] <- 'invented'; fails('classification', df, message='Invalid convergence_type')
fails('duplicate', rbind(base, base[1,]), message='Duplicate')
df <- base; df$side[1] <- 'none'; df$asr_path_score[1] <- 0
fails('conflicting_none', df, message='Conflicting none-side')
df <- base; df$asr_path_score[3] <- .1; fails('positive_none', df, message='Conflicting none-side')
rec <- record_for(base); rec$complete <- FALSE; fails('incomplete', record=rec, message='Invalid observed')
rec <- record_for(base); rec$batches[[1]]$score_rule <- 'old'; fails('rule', record=rec, message='incompatible')
rec <- record_for(base); rec$batches[[1]]$posterior_threshold <- 2; fails('threshold', record=rec, message='Invalid recorded')
rec <- record_for(base); rec$genes$`001.2`$sha256 <- 'wrong'; fails('ledger_coverage', record=rec, message='coverage differs')
fails('truncated_master', base[-1:-2,], record=record_for(base), message='gene coverage differs')
df <- base; df$msa_pos[3] <- 99; fails('changed_position', df, record=record_for(base), message='position coverage differs')
rec <- record_for(base); second <- rec$batches[[1]]; second$genes <- structure(list(), names=character())
second$posterior_threshold <- .3; rec$batches[[2]] <- second
fails('inconsistent_batches', record=rec, message='Inconsistent scientific')
rec <- record_for(base); rec$batches[[2]] <- rec$batches[[1]]
fails('overlapping_batches', record=rec, message='Overlapping genes')
rec <- record_for(base); rec$excluded_genes <- list('G.3'); fails('excluded_overlap', record=rec, message='both scored and excluded')
d <- write_case('missing_stamp'); unlink(file.path(d,'observed_scoring.json'))
stopifnot(inherits(tryCatch(load_asr_diagnostics(d),error=identity),'error'))

# Five full renders cover empty data, absent optional columns, unavailable agreement,
# all categories including diagnostic zeros, and an optional current threshold.
for (case in c('mixed', 'empty', 'required_only', 'no_agreement', 'recorded_only')) {
  df <- base
  if (case == 'empty') df <- df[FALSE,]
  if (case == 'required_only') df <- df[c('gene','msa_pos','side','asr_path_score')]
  if (case == 'no_agreement') df$derived_agreement <- NA_real_
  d <- write_case(paste0('render_',case), df)
  file.copy(list.files(file.path(root,'subworkflows/ASR_DIAGNOSTICS/local'), full.names=TRUE), d, recursive=TRUE)
  file.copy(file.path(root,'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'), file.path(d,'src'))
  requested <- if (case == 'recorded_only') NULL else .1
  e <- new.env(parent=globalenv())
  rmarkdown::render(file.path(d,'9.ASR_diagnostics.Rmd'),
    params=list(disambig_dir=d, output_dir=d, posterior_threshold=requested), quiet=TRUE, envir=e)
  html <- paste(readLines(file.path(d,'9.ASR_diagnostics.html'),warn=FALSE), collapse='\n')
  stopifnot(grepl('Helvetica Neue',html,fixed=TRUE),grepl('Recorded posterior threshold',html,fixed=TRUE),
            e$posterior_threshold == .2, e$threshold_mismatch == !is.null(requested),
            length(list.files(file.path(d,'tsv'),pattern='tsv$')) == 2,
            !any(is.infinite(e$summary_table$Mean)))
  if (case == 'empty') stopifnot(nrow(e$summary_table)==0, !length(list.files(file.path(d,'plots'))))
  if (case %in% c('required_only','no_agreement')) stopifnot(!file.exists(file.path(d,'plots/derived_agreement_distribution.png')))
  if (case == 'required_only') stopifnot(all(is.na(e$agree_tab$`Convergent (%)`)))
  if (case == 'mixed') {
    none <- e$agree_tab[e$agree_tab$Side == 'none',]
    stopifnot(abs(sum(none[c('No change (%)','Symplesiomorphic (%)','Cross hypothesis (%)')])-100)<1e-10,
              identical(e$p_score$layers[[1]]$aes_params$fill, e$.tableau10[1]))
  }
  cat('PASS full render:', case, '\n')
}
cat('ASR diagnostics contract and render checks complete.\n')
