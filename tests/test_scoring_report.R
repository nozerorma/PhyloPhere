# Rscript tests/test_scoring_report.R <repo root> <scratch directory>
# Render the whole report so genomic setup, plots, tables and exports run together.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
scratch <- normalizePath(args[2])

genes <- data.frame(
  Gene = c("g1", "g2", "g3"), n_positions = c(2L, 1L, 0L),
  gene_caas_score = c(.8, .4, 0), gene_caas_score_top = c(.8, .4, 0),
  gene_caas_score_bottom = c(.2, 0, 0), caas_p.gen = c(.01, .5, NA_real_)
)
positions <- data.frame(
  Gene = c("g1", "g1", "g2", "unmapped"), Position = c(0L, 1L, 0L, 0L),
  side = c("top", "bottom", "top", "top"), CAAS_score = c(.8, .2, .4, .9),
  p.emp = c(.01, .5, .2, .1), p.adj_bh = c(.04, .5, .4, .2)
)
coords <- data.frame(
  Gene = c("g3", "g2", "g1"), chr = c("chrX", "chr2", "chr1"),
  start = c(500L, 200L, 100L), end = c(600L, 300L, 200L), length = 101L
)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)

for (scenario in c("mapped", "unmapped_positions", "no_coordinates")) {
  wd <- file.path(scratch, scenario)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(root, "subworkflows/SCORING/local/11.Scoring_report.Rmd"), wd, overwrite = TRUE)
  dir.create(file.path(wd,"src"),showWarnings=FALSE)
  file.copy(file.path(root,"subworkflows/SCORING/local/src/report_evidence.R"),file.path(wd,"src"),overwrite=TRUE)
  jsonlite::write_json(list(gene_stat="max",nominal_threshold=.05,fdr_threshold=.05), file.path(wd,"gene_evidence_metadata.json"),auto_unbox=TRUE)
  write_tsv(genes, file.path(wd, "gene_scores.tsv"))
  pos <- positions
  pos$p.emp_fact <- .0000001  # legacy input must not leak into reports or exports
  pos$p.adj_bh_fact <- .000001
  if (scenario == "unmapped_positions") pos$Gene <- paste0("missing", seq_len(nrow(pos)))
  write_tsv(pos, file.path(wd, "position_scores.tsv"))
  write_tsv(data.frame(Gene = character(), correlation = numeric()), file.path(wd, "gene_correlations.tsv"))
  write_tsv(coords, file.path(wd, "coordinates.tsv"))
  env <- new.env(parent = globalenv())
  rmarkdown::render(file.path(wd, "11.Scoring_report.Rmd"), quiet = TRUE, envir = env,
    params = list(position_scores_file = "position_scores.tsv", gene_scores_file = "gene_scores.tsv",
      gene_corr_file = "gene_correlations.tsv", output_dir = wd, window_size_bp = 100000,
      genomic_info_file = if (scenario == "no_coordinates") NULL else "coordinates.tsv"),
    output_file = "report.html")
  stopifnot(file.exists(file.path(wd, "report.html")),
            file.exists(file.path(wd, "position_scores_annotated.tsv")))
  exported <- read.delim(file.path(wd, "position_scores_annotated.tsv"))
  stopifnot(!any(grepl("fact", names(exported))),
    !any(grepl("p[.]emp_fact|p[.]adj_bh_fact", readLines(file.path(wd, "report.html")))))
  if (scenario != "no_coordinates") {
    wins <- read.delim(file.path(wd, "genomic_window_accum_gene_caas.tsv"))
    stopifnot(sum(wins$n_genes) == 3L,
      isTRUE(all.equal(sort(wins$mean_score),c(0,.4,.8))),
      !file.exists(file.path(wd,"genomic_window_accum_positions.tsv")))
  } else stopifnot(!exists("gene_genomic", envir=env, inherits=FALSE))
  cat("PASS scoring report:", scenario, "\n")
}
