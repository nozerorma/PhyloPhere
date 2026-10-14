# Streaming readers of perm_pos_cycle_caas.tsv.gz (scoring_compute.R) give what the in-memory read gave.
# Run: Rscript validation/unification/test_null_streaming.R   (from the repository root)
suppressPackageStartupMessages({library(readr); library(dplyr)})
src <- readLines("subworkflows/SCORING/local/src/scoring_compute.R")
a <- grep("^  \\.null_header <- ", src); b <- grep("^  \\.cycles_in_file <- \\.meta_val", src)
stream_code <- paste(src[a:b], collapse = "\n")

set.seed(11)
tmp <- tempfile("nullstream_", tmpdir = file.path(Sys.getenv("HOME"), ".cc")); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
make_null <- function(path, n_genes = 120, n_cycles = 40, agg = "mean_species", interleave_gene = FALSE) {
  rows <- vector("list", n_genes)
  for (g in seq_len(n_genes)) {
    n <- sample(5:600, 1)
    d <- data.frame(Gene = sprintf("G%03d", g), Position = sample(1:900, n, TRUE),
                    side = sample(c("top", "bottom", "none"), n, TRUE, prob = c(.4, .4, .2)),
                    cycle = sprintf("b_%d", sample(1:n_cycles, n, TRUE)),
                    caas_score = sprintf("%.17g", runif(n)),
                    score_aggregation = agg, stringsAsFactors = FALSE)
    d$caas_score[sample.int(n, floor(n * .6))] <- ""          # detections without a score
    d$caas_score[sample.int(n, floor(n * .05))] <- "0"        # a zero is not positive
    d <- rbind(d, within(d[sample.int(n, floor(n / 3)), ], side <- ifelse(side == "top", "bottom", "top")))  # both sides
    rows[[g]] <- d
  }
  x <- bind_rows(rows)
  if (interleave_gene) x <- x[sample(nrow(x)), ]
  write_tsv(x, path, na = "")
}

reference <- function(f, score_aggregation) {
  cyc_caas <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = col_guess(), caas_score = col_character()))
  list(pooled = cyc_caas %>% mutate(Position = as.integer(Position), caas_score = suppressWarnings(as.numeric(caas_score))) %>%
         filter(!is.na(caas_score) & caas_score > 0) %>% group_by(Gene, Position, cycle) %>%
         summarise(caas_max = max(caas_score), .groups = "drop"),
       cycles = unique(as.character(cyc_caas$cycle)), agg = unique(as.character(cyc_caas$score_aggregation)))
}
streamed <- function(f, score_aggregation) {
  caas_pos_cycle_caas_file <- f
  eval(parse(text = stream_code))
  list(pooled = cyc_pooled, cycles = .cycles_in_file, agg = .null_agg)
}

f <- file.path(tmp, "null.tsv.gz"); make_null(f)
ref <- reference(f, "mean_species"); got <- streamed(f, "mean_species")
key <- function(x) x %>% arrange(Gene, Position, cycle) %>% as.data.frame()
stopifnot(identical(key(ref$pooled)$caas_max, key(got$pooled)$caas_max),
          identical(key(ref$pooled)[, 1:3], key(got$pooled)[, 1:3]),
          setequal(ref$cycles, got$cycles), identical(ref$agg, got$agg))
cat(sprintf("pooled null: %d (Gene, Position, cycle) rows identical, %d cycles, rule '%s'\n", nrow(got$pooled), length(got$cycles), got$agg))

# A null of another score rule is refused as before
f2 <- file.path(tmp, "other.tsv.gz"); make_null(f2, n_genes = 10, agg = "max")
r <- tryCatch(streamed(f2, "mean_species"), error = function(e) conditionMessage(e)); stopifnot(is.character(r), grepl("scored with the rule", r))
# Rows of a gene that are not together are refused
f3 <- file.path(tmp, "mixed.tsv.gz"); make_null(f3, n_genes = 10, interleave_gene = TRUE)
r <- tryCatch(streamed(f3, "mean_species"), error = function(e) conditionMessage(e)); stopifnot(is.character(r), grepl("not contiguous", r))
cat("rule mismatch and non-contiguous genes refused\nOK\n")
