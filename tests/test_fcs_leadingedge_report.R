# Rscript tests/test_fcs_leadingedge_report.R <repo root> <scratch directory>
# Uses precomputed enrichment, as MODULE_FCS_RER does, and renders the full report.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
scratch <- normalizePath(args[2])
local_dir <- file.path(root, "subworkflows/ENRICHMENT/local")
genes <- paste0("g", 1:4)
stats <- data.frame(gene = genes, score_global = c(4, 3, -3, -4),
                    score_top = c(4, 3, 0, 0), score_bottom = c(0, 0, -3, -4),
                    score_custom = c(4, 3, -3, -4))
annot <- data.frame(gene = genes,
  caas_nominal = c(TRUE, TRUE, FALSE, FALSE), caas_adjusted = c(TRUE, FALSE, FALSE, FALSE),
  flag_rer_acc = c(TRUE, TRUE, FALSE, FALSE), flag_rer_decc = c(FALSE, FALSE, TRUE, TRUE),
  flag_rer_acc_adjusted = c(TRUE, FALSE, FALSE, FALSE),
  flag_rer_decc_adjusted = c(FALSE, FALSE, TRUE, FALSE),
  flag_custom = c(TRUE, FALSE, FALSE, FALSE))
enrich <- data.frame(ranking = c("global", "top", "bottom", "custom"),
  database = "test", pathway = "setA", stat = c(.4, .4, -.4, .4),
  pval = .001, p.adj = .002, p.perm = .01, q.perm = .02, num.genes = 4,
  gene.vals = "g1:4,g2:3,g3:2,g4:1", sig_perm = TRUE,
  direction = c("positive", "positive", "negative", "positive"), redlich_level = .05)

for (scenario in c("significant", "no_significant_sets")) {
  wd <- file.path(scratch, scenario)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(local_dir, "12.FCS_general_report.Rmd"), wd, overwrite = TRUE)
  dir.create(file.path(wd, "src"), showWarnings = FALSE)
  file.copy(file.path(local_dir, "src", c("fcs_enrich.R", "fcs_cluster.R", "evidence_flags.R")),
            file.path(wd, "src"), overwrite = TRUE)
  dir.create(file.path(wd, "gmts"), showWarnings = FALSE)
  writeLines(c(paste(c("setA", "Test pathway", genes), collapse = "\t"),
               paste(c("setB", "Other pathway", genes[1:2]), collapse = "\t")),
             file.path(wd, "gmts/test.gmt"))
  writeLines(genes, file.path(wd, "background.txt"))
  readr::write_tsv(stats, file.path(wd, "stats.tsv"))
  readr::write_tsv(annot, file.path(wd, "annot.tsv"))
  ev <- enrich
  if (scenario == "no_significant_sets") {
    ev$q.perm <- 1
    ev$sig_perm <- FALSE
    ev$redlich_level <- NA_real_
  }
  readr::write_tsv(ev, file.path(wd, "enrich.tsv"))
  env <- new.env(parent = globalenv())
  rmarkdown::render(file.path(wd, "12.FCS_general_report.Rmd"), quiet = TRUE, envir = env,
    params = list(stats_file = "stats.tsv", annot_file = "annot.tsv", universe_file = "background.txt",
      gmt_dir = file.path(wd, "gmts"), enrich_file = "enrich.tsv", score_source = "rer", num_g = 1),
    output_file = "report.html")
  stopifnot(file.exists(file.path(wd, "report.html")))
  if (scenario == "significant") {
    comp <- readr::read_tsv(file.path(wd, "fcs_results/fcs_leading_edge_composition.tsv"),
                            show_col_types = FALSE, name_repair = "minimal")
    stopifnot(!anyNA(names(comp)), all(nzchar(names(comp))), !anyDuplicated(names(comp)),
      setequal(comp$ranking, c("global", "top", "bottom", "custom")),
      all(comp$n_le == 2), comp$`rer_acc BH`[comp$ranking == "global"] == 50,
      comp$`rer_acc BH`[comp$ranking == "top"] == 50,
      comp$`rer_decc BH`[comp$ranking == "bottom"] == 50,
      comp$flag_custom[comp$ranking == "custom"] == 50)
  } else {
    stopifnot(nrow(env$le_comp_export) == 0L)
  }
  cat("PASS FCS leading-edge report:", scenario, "\n")
}
