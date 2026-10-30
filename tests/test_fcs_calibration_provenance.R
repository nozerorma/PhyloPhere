# The calibrated FCS refuses a universe other than the one of the calibration and takes its gene batches from the
# universe as it is read (unique, trimmed, non-empty), as the calibration does. Run from the repository root.
suppressPackageStartupMessages({library(dplyr); library(readr)})
src <- file.path(normalizePath("."), "subworkflows/ENRICHMENT/local/src")
for (f in c("gene_p.R", "caas_blocks.R")) source(file.path(src, f))
wd <- tempfile("fcs_prov_"); dir.create(wd); setwd(wd)
cli <- function(...) {
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c(file.path(src, "fcs_blocks.R"), ...), stdout = TRUE, stderr = TRUE))
  list(ok = is.null(attr(out, "status")), out = paste(out, collapse = "\n"))
}
set.seed(3)
genes <- paste0("g", 1:12)
writeLines(c(genes[1:6], "", genes[1:6], genes[7:12]), "universe.txt")      # blank and repeated lines: 12 genes as read
cycles <- paste0("c", 1:20)
raw <- matrix(sample(0:2, 12 * 20, replace = TRUE), 12, dimnames = list(genes, cycles))
observed <- setNames(sample(0:3, 12, replace = TRUE), genes)
gp <- gene_p(raw, observed)
G <- 5; batches <- caas_batches(genes, G)
for (k in seq_along(batches)) {
  rows <- batches[[k]]
  write_tsv(tibble(gene = rows, nlog10_obs = gp$nlog10_obs[rows]), sprintf("obs_global_%04d.tsv", k))
  saveRDS(list(rows = rows, cols = cycles, values = gp$nlog10_null[rows, , drop = FALSE]), caas_slice_path(".", "global_gp", k, "blk1"))
}
md5 <- unname(tools::md5sum("universe.txt"))
jsonlite::write_json(list(universe_md5 = md5, genes = 12L, available = TRUE), "calibration_metadata.json", auto_unbox = TRUE)
dir.create("gmts"); writeLines(c(paste(c("A", "a", genes[1:6]), collapse = "\t"), paste(c("B", "b", genes[5:10]), collapse = "\t"),
                                 paste(c("C", "c", genes[7:12]), collapse = "\t")), "gmts/t.gmt")
observed_step <- function(universe = "universe.txt", meta = "calibration_metadata.json")
  cli("observed", "--obs-dir", ".", "--universe", universe, "--rankings", "global", "--gmt-dir", "gmts", "--num-g", 3,
      "--max-g", 0, "--metadata", meta, "--out", "observed.rds")
ok <- observed_step(); stopifnot(ok$ok)
# another universe (reversed order, so another md5) is rejected with the reason
writeLines(rev(genes), "universe_other.txt")
bad <- observed_step("universe_other.txt"); stopifnot(!bad$ok, grepl("universe_md5", bad$out))
# the gene count of the calibration must match as well
jsonlite::write_json(list(universe_md5 = md5, genes = 11L, available = TRUE), "meta_count.json", auto_unbox = TRUE)
bad <- observed_step(meta = "meta_count.json"); stopifnot(!bad$ok, grepl("gene count", bad$out))
# observed gene-p batches must cover the universe
file.rename("obs_global_0003.tsv", "obs_global_0003.bak")
bad <- observed_step(); stopifnot(!bad$ok, grepl("cover the FCS universe", bad$out))
file.rename("obs_global_0003.bak", "obs_global_0003.tsv")
# null_block derives its batches from the universe as read: 3 batches of 5 genes, not the 13 non-empty lines of the file
nb <- cli("null_block", "--observed", "observed.rds", "--slices", ".", "--universe", "universe.txt", "--block", "blk1",
          "--batch-size", G, "--rank", "global", "--alt", "greater", "--gmt-dir", "gmts", "--num-g", 3, "--max-g", 0, "--out", "counts.tsv")
stopifnot(nb$ok)
counts <- read_tsv("counts.tsv", show_col_types = FALSE)
stopifnot(nrow(counts) > 0, all(counts$valid == 20))
cat("PASS calibrated FCS: universe provenance and gene batches\n")
