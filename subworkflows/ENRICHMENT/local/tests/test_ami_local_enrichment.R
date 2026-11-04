# The memoized local enrichment of the AMI report returns exactly what the per-call computation returns,
# for several modules, backgrounds and the no-background fallback, and does not mix cached backgrounds.
args <- commandArgs(TRUE); root <- args[1]
suppressPackageStartupMessages(library(dplyr))
src <- readLines(file.path(root, 'subworkflows/ENRICHMENT/local/12.AMI_analysis.Rmd'))
from <- grep('^\\.enrichment_memo <- new.env', src); to <- grep('^# ── Palette and file helpers', src)
stopifnot(length(from) == 1L, length(to) == 1L, from < to)

# Reference: the computation that ran once per module.
reference <- function(fg_ids, bg_ids = NULL, fdr_threshold = .description_fdr) {
  if (length(fg_ids) == 0) return(tibble::tibble())
  if (is.null(bg_ids) || length(bg_ids) < 5) bg_ids <- unique(enrich_terms_db$string_protein_id)
  fg_ids <- intersect(fg_ids, bg_ids)
  if (length(fg_ids) == 0) return(tibble::tibble())
  bg_data <- enrich_terms_db %>% dplyr::filter(string_protein_id %in% bg_ids)
  term_bg_counts <- bg_data %>% dplyr::group_by(category, term, description) %>%
    dplyr::summarise(M = dplyr::n(), .groups = "drop") %>% dplyr::filter(M >= 2)
  fg_data <- enrich_terms_db %>% dplyr::filter(string_protein_id %in% fg_ids)
  fg_counts <- fg_data %>% dplyr::group_by(category, term, description) %>% dplyr::summarise(k = dplyr::n(), .groups = "drop")
  N <- length(bg_ids); n <- length(fg_ids)
  res <- term_bg_counts %>% dplyr::inner_join(fg_counts, by = c("category", "term", "description")) %>%
    dplyr::mutate(pvalue = phyper(k - 1, M, N - M, n, lower.tail = FALSE)) %>%
    dplyr::group_by(category) %>% dplyr::mutate(fdr = p.adjust(pvalue, method = "fdr")) %>% dplyr::ungroup() %>%
    dplyr::filter(fdr <= fdr_threshold) %>% dplyr::arrange(category, fdr)
  if (nrow(res) == 0) return(tibble::tibble())
  fg_mapped <- fg_data %>% dplyr::mutate(gene = sapply(string_protein_id, function(id) { s <- global_id2gene[[id]]; if (is.null(s)) id else s }))
  term_genes <- fg_mapped %>% dplyr::group_by(category, term) %>%
    dplyr::summarise(preferredNames = paste(sort(unique(gene)), collapse = ","), .groups = "drop")
  res %>% dplyr::left_join(term_genes, by = c("category", "term")) %>%
    dplyr::rename(number_of_genes = k, number_of_genes_in_background = M, p_value = pvalue)
}

set.seed(7)
proteins <- sprintf("9606.ENSP%05d", 1:3000)
n_terms <- 400
member <- lapply(seq_len(n_terms), function(i) sample(proteins, sample(5:150, 1)))
enrich_terms_db <- tibble::tibble(
  string_protein_id = unlist(member),
  category = rep(c("Process", "Function", "KEGG")[(seq_len(n_terms) %% 3) + 1], lengths(member)),
  term = rep(sprintf("T%03d", seq_len(n_terms)), lengths(member)),
  description = rep(sprintf("term %03d", seq_len(n_terms)), lengths(member)))
global_id2gene <- new.env(); for (p in proteins[1:500]) global_id2gene[[p]] <- sub("9606.ENSP", "SYM", p)
.description_fdr <- 0.10
eval(parse(text = src[from:(to - 1)]))   # defines memoized(), database_proteins(), background_term_counts(), get_enrichment_local()

bg_a <- sample(proteins, 2000); bg_b <- sample(proteins, 1200)
n_nonempty <- 0L
for (bg in list(bg_a, bg_b, NULL, bg_a, c(bg_a, bg_a[1:10]))) {      # the repeats hit the cache
  for (k in 1:6) {
    m <- member[[sample(n_terms, 1)]]
    fg <- if (k %% 2) sample(m, min(20, length(m))) else sample(proteins, 40)
    got <- get_enrichment_local(fg, bg); want <- reference(fg, bg)
    stopifnot(identical(got, want))
    n_nonempty <- n_nonempty + (nrow(want) > 0)
  }
}
stopifnot(n_nonempty >= 5L)                                          # the comparison is not vacuous
stopifnot(identical(get_enrichment_local(character(0), bg_a), tibble::tibble()))
stopifnot(sum(startsWith(ls(.enrichment_memo), "terms_")) == 3L)     # bg_a, bg_b, all proteins; duplicates share a key
cat('AMI_LOCAL_ENRICHMENT_PASS\n')
