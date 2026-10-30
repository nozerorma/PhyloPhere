# Run from the repository root: Rscript tests/test_fcs_rer_missing.R
# Optional first argument: saved marine_complete results directory.
repo <- normalizePath(".")
src <- file.path(repo, "subworkflows/ENRICHMENT/local/src")
source(file.path(src, "evidence_flags.R"))
source(file.path(src, "fcs_enrich.R"))

local({
  tmp <- tempfile("fcs_rer_missing_")
  dir.create(tmp)
  old <- setwd(tmp)
  on.exit({ setwd(old); unlink(tmp, recursive = TRUE) })
  same <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE)))
  key <- function(x) paste(x$ranking, x$database, x$pathway, sep = "|")
  compare <- function(a, b) {
    stopifnot(setequal(key(a), key(b)))
    b <- b[match(key(a), key(b)), , drop = FALSE]
    stopifnot(!"sig_top" %in% names(a), !"sig_top" %in% names(b))
    for (n in c("stat", "pval", "p.adj", "p.perm", "q.perm", "num.genes", "sig_perm", "redlich_level")) same(as.numeric(a[[n]]), as.numeric(b[[n]]))
  }

  genes <- sprintf("g%02d", 1:40)
  res <- data.frame(Rho = rep(c(-.5, .5), 20), P = seq(.001, .9, length.out = 40), row.names = genes)
  res$Rho[1] <- NA_real_; res$P[2] <- NA_real_; res$Rho[3] <- 0
  scores <- setNames(sign(res$Rho) * -log10(res$P), genes)
  universe <- c(genes, "background_only")
  rer <- fcs_build_vals(scores, universe, score_source = "rer")
  same(rer, RERconverge::getStat(res))
  stopifnot(identical(names(rer), genes[-c(1, 2)]), rer[[genes[3]]] == 0,
            length(fcs_build_vals(setNames(c(NA_real_, NaN), genes[1:2]), universe, score_source = "rer")) == 0)
  caas <- fcs_build_vals(scores, universe)
  same(caas, fcs_build_vals(scores, universe, score_source = "caas"))
  stopifnot(identical(names(caas), universe), all(caas[c(genes[1:3], "background_only")] == 0))

  dir.create("gmts")
  writeLines(c(paste(c("A", "set A", genes[1:15]), collapse = "\t"),
               paste(c("B", "set B", genes[12:28]), collapse = "\t"),
               paste(c("C", "set C", genes[25:40]), collapse = "\t")), "gmts/test.gmt")
  invisible(capture.output(suppressMessages(gmts <- fcs_load_gmts("gmts"))))
  set.seed(42)
  null <- matrix(rnorm(40 * 30), 40, dimnames = list(genes, NULL))
  null[1, 1:5] <- NA_real_; null[2, 6:10] <- NA_real_; null[3, 11:15] <- 0
  saveRDS(list(corStat = null), "perms.rds")
  readr::write_tsv(data.frame(gene = genes, score_global = scores), "stats.tsv")
  writeLines(universe, "universe.txt")
  reference <- fcs_run_all(list(global = RERconverge::getStat(res)), gmts, perms_file = "perms.rds")
  stopifnot(nrow(reference) == 3, !anyNA(reference$p.perm))
  # Execute the report's actual fallback expression, independently of HTML rendering.
  rmd <- readLines(file.path(src, "../11.FCS_general_report.Rmd"))
  # The if/else of the run chunk: the last block that starts on enrich_file, up to the sourcing of fcs_selection.R.
  start <- tail(grep('^if \\(!is.null\\(params\\$enrich_file\\)', rmd), 1)
  end <- grep('^source\\(file.path\\(dirname\\(src\\), "fcs_selection.R"\\)\\)', rmd)
  fallback <- parse(text = rmd[start:(end - 1)])
  stats <- readr::read_tsv("stats.tsv", show_col_types = FALSE)
  score_cols <- "score_global"
  params <- list(enrich_file = NULL, score_source = "rer", num_g = 10, max_g = 0,
                 perms_file = "perms.rds", fdr_thr = .15, pperm_thr = .025, analytic_levels = "0.05,0.1,0.15", gene_stat = "max")
  eval(fallback)
  compare(enrich, reference)
  # A signed test remains bilateral when this sample happens to have only positive scores.
  positive <- rer - min(rer) + 1
  stopifnot(identical(fcs_alternative(positive, "rer"), "two.sided"),
            identical(fcs_alternative(positive, "caas"), "greater"))
  signed_ref <- fcs_run_ranking(positive, gmts, alternative = "two.sided")
  forced <- fcs_observed(list(global = positive), gmts, score_source = "rer")
  stopifnot(identical(forced$alts$global, "two.sided"))
  same(forced$enrich$stat, signed_ref$stat)
  same(forced$enrich$pval, signed_ref$pval)
  # The report fallback also carries the rer source through to the test.
  stats$score_global <- positive[stats$gene]
  readr::write_tsv(stats, "stats.tsv")
  expected <- fcs_run_all(list(global = positive), gmts, perms_file = "perms.rds", score_source = "rer")
  eval(fallback)
  compare(enrich, expected)
  stopifnot(anyNA(null), any(null == 0, na.rm = TRUE))
  cat("PASS: native getStat membership, measured zeros, all-missing scores, CAAS unchanged and report fallback\n")

  marine <- commandArgs(TRUE)
  if (length(marine)) {
    run <- normalizePath(marine[1])
    observed <- readRDS(file.path(run, "rerconverge/rer_results/marine.binary.output"))
    perms <- file.path(run, "rerconverge/rer_results/marine.binary.perms.rds")
    invisible(capture.output(suppressMessages(gmts <- fcs_load_gmts(file.path(run, "core_inputs/gmt")))))
    scores <- setNames(sign(observed$Rho) * -log10(pmax(observed$P, 1e-300)), rownames(observed))
    before <- fcs_run_all(list(global = fcs_build_vals(scores, rownames(observed))), gmts, perms_file = perms)
    after <- fcs_run_all(list(global = fcs_build_vals(scores, rownames(observed), score_source = "rer")), gmts, perms_file = perms)
    native <- fcs_run_all(list(global = RERconverge::getStat(observed)), gmts, perms_file = perms)
    compare(after, native)
    stopifnot(nrow(after) == 133, sum(before$sig_perm) == 0, sum(after$sig_perm) == 1)
    cat(sprintf("PASS: marine, %d observed genes -> %d scored genes, %d sets, sig_perm %d -> %d\n",
                length(scores), length(RERconverge::getStat(observed)), nrow(after),
                sum(before$sig_perm), sum(after$sig_perm)))
  }
})
