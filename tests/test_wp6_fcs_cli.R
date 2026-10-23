# WP6 steps as separate Rscript calls (fcs_blocks.R) equal fcs_run_all on the whole null with the -log10 p of GENE_P.
# Run from repository root; compares physical CLI steps with the whole-memory AUC contract.
new_dir <- normalizePath('subworkflows/ENRICHMENT/local/src')
old_file <- file.path(new_dir, 'fcs_enrich.R')
suppressPackageStartupMessages({ library(RERconverge); library(readr) })
source(file.path(new_dir, "evidence_flags.R"))
env <- new.env(); sys.source(file.path(new_dir, "fcs_enrich.R"), envir = env)
sys.source(file.path(new_dir, "gene_p.R"), envir = env); sys.source(file.path(new_dir, "caas_blocks.R"), envir = env)
env_old <- new.env(); sys.source(old_file, envir = env_old)
cli <- file.path(new_dir, "fcs_blocks.R")
run_step <- function(...) { out <- system2(Sys.getenv("RSCRIPT", "Rscript"), c(cli, ...), stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status"); if (!is.null(st) && st != 0) stop(paste(out, collapse = "\n")); invisible(out) }

set.seed(3)
universe <- paste0("G", sprintf("%03d", 1:200))
sets <- lapply(1:40, function(i) sample(universe, sample(10:40, 1)))
root <- tempfile("cli"); dir.create(root)
gmt_dir <- file.path(root, "gmts"); dir.create(gmt_dir)
writeLines(vapply(seq_along(sets), function(i) paste(c(paste0("S", i), "desc", sets[[i]]), collapse = "\t"), ""),
           file.path(gmt_dir, "dbA.gmt"))
writeLines(universe, file.path(root, "universe.txt"))
cycles_by_block <- list(blk_0001 = paste0("b_", 1:4), blk_0002 = paste0("b_", 5:9), blk_0003 = paste0("b_", 10:11))
all_cycles <- unlist(cycles_by_block)
raw <- matrix(abs(rnorm(200 * length(all_cycles))), 200, length(all_cycles), dimnames = list(universe, all_cycles))
raw[sample(length(raw), 800)] <- 0
obs_stat <- setNames(abs(rnorm(200)) * 0.5 + 0.2, universe); obs_stat[1:10] <- 0
block_mats <- lapply(cycles_by_block, function(cy) { m <- raw[, cy, drop = FALSE]; m[sample(universe, 20), ] <- 0; m })
full_ref <- raw
for (b in names(cycles_by_block)) full_ref[, cycles_by_block[[b]]] <- block_mats[[b]]
ref_gp <- env$gene_p(full_ref, obs_stat)
rankings <- list(global = env$fcs_build_vals(setNames(ref_gp$nlog10_obs, universe), universe))
perms_file <- file.path(root, "perms.rds")
saveRDS(list(caas_corStat_byrank = list(global = ref_gp$nlog10_null), gene_stat = "size_adj_max"), perms_file)
gmts_ref <- env_old$fcs_load_gmts(gmt_dir)
ref <- env_old$fcs_run_all(rankings, gmts_ref, num_g = 10, max_g = 0, perms_file = perms_file, fdr_thr = 0.15,
                           p_perm_thr = 0.025, gene_stat = "size_adj_max")

ok <- TRUE
for (G in c(200)) {
  wd <- file.path(root, sprintf("g%d", G)); dir.create(wd); slices <- file.path(wd, "slices"); dir.create(slices)
  uni <- file.path(root, "universe.txt")
  for (b in names(block_mats)) {
    mf <- file.path(wd, paste0(b, "_matrices.rds"))
    saveRDS(list(genes = universe, cycles = cycles_by_block[[b]], gene_stat = "size_adj_max",
                 stat = list(global = block_mats[[b]])), mf)
    run_step("transpose", "--matrices", mf, "--universe", uni, "--batch", G, "--rank", "global", "--block", b, "--out", slices)
  }
  obs_tab <- file.path(wd, "obs.tsv")
  write.table(data.frame(gene = universe, gene_caas_score = obs_stat), obs_tab, sep = "\t", quote = FALSE, row.names = FALSE)
  nb <- ceiling(200 / G)
  for (k in seq_len(nb)) run_step("gene_p", "--batch", k, "--rank", "global", "--blocks", paste(names(block_mats), collapse = ","),
                                  "--slices", slices, "--obs", obs_tab, "--obs-col", "gene_caas_score", "--out", wd)
  # diagnostics: each stage against the in-memory reference
  og <- do.call(rbind, lapply(list.files(wd, pattern = "^obs_global_", full.names = TRUE), read.delim))
  cat(sprintf("  stage gene_p: observed -log10 p max |diff| %.3g over %d genes\n",
              max(abs(og$nlog10_obs - ref_gp$nlog10_obs[og$gene])), nrow(og)))
  ref_obs <- env$fcs_observed(rankings, gmts_ref, num_g = 10, max_g = 0)
  cat(sprintf("  stage observed (in memory, previous rows): %d rows\n", nrow(ref_obs$enrich)))
  run_step("observed", "--obs-dir", wd, "--universe", uni, "--rankings", "global", "--gmt-dir", gmt_dir, "--out", file.path(wd, "observed.rds"))
  for (b in names(block_mats)) run_step("null_block", "--observed", file.path(wd, "observed.rds"), "--slices", slices, "--universe", uni,
      "--block", b, "--batches", paste(seq_len(nb), collapse = ","), "--rank", "global", "--alt", "greater", "--gmt-dir", gmt_dir,
      "--out", file.path(wd, paste0("counts_", b, ".tsv")))
  run_step("close", "--observed", file.path(wd, "observed.rds"),
           "--counts", paste(file.path(wd, paste0("counts_", names(block_mats), ".tsv")), collapse = ","), "--out", file.path(wd, "enrich.tsv"))
  cl <- do.call(rbind, lapply(names(block_mats), function(b) readr::read_tsv(file.path(wd, paste0("counts_", b, ".tsv")), show_col_types = FALSE)))
  cc <- dplyr::group_by(cl, database, pathway) %>% dplyr::summarise(count = sum(count), valid = sum(valid), .groups = "drop")
  mem_mats <- lapply(names(cycles_by_block), function(b) ref_gp$nlog10_null[, cycles_by_block[[b]], drop = FALSE])
  mem <- env$fcs_null_counts_by_matrices("global", ref_obs$enrich[ref_obs$enrich$ranking == "global", ], mem_mats, gmts_ref, "greater", num_g = 10, max_g = 0)
  mm <- match(paste(mem$database, mem$pathway), paste(cc$database, cc$pathway))
  cat(sprintf("  stage counts vs in memory: count max |diff| %.3g, valid max |diff| %.3g\n",
              max(abs(cc$count[mm] - mem$count)), max(abs(cc$valid[mm] - mem$valid))))
  got <- readr::read_tsv(file.path(wd, "enrich.tsv"), show_col_types = FALSE)
  ob <- readRDS(file.path(wd, "observed.rds"))
  mo <- match(paste(ref_obs$enrich$ranking, ref_obs$enrich$database, ref_obs$enrich$pathway), paste(ob$ranking, ob$database, ob$pathway))
  cat(sprintf("  stage observed: rows %d vs %d, stat max |diff| %.3g, pval max |diff| %.3g\n", nrow(ob), nrow(ref_obs$enrich),
              max(abs(ob$stat[mo] - ref_obs$enrich$stat)), max(abs(ob$pval[mo] - ref_obs$enrich$pval))))
  cat(sprintf("  stage counts: %s\n", paste(sprintf("%s=%d", names(block_mats), vapply(file.path(wd, paste0("counts_", names(block_mats), ".tsv")), function(f) nrow(readr::read_tsv(f, show_col_types = FALSE)), 1L)), collapse = " ")))
  key <- function(d) paste(d$ranking, d$database, d$pathway)
  m <- match(key(ref), key(got))
  same_keys <- !anyNA(m) && nrow(got) == nrow(ref)
  pdiff <- max(abs(got$p.perm[m] - ref$p.perm), na.rm = TRUE)
  same_na <- identical(is.na(got$p.perm[m]), is.na(ref$p.perm))
  sig_same <- identical(as.logical(got$sig_perm[m]), as.logical(ref$sig_perm))
  stopifnot(all(cc$count[mm] == mem$count), all(cc$valid[mm] == mem$valid))
  pass <- same_keys && pdiff <= 1e-12 && same_na && sig_same && nrow(ref) > 0
  cat(sprintf("batch %3d genes, 3 blocks, CLI steps: sets %d, p.perm max |diff| %.3g, NA equal %s, significance equal %s -> %s\n",
              G, nrow(ref), pdiff, same_na, sig_same, if (pass) "PASS" else "FAIL"))
  ok <- ok && pass
}
cat(if (ok) "ALL PASS\n" else "FAILURES\n"); quit(status = if (ok) 0 else 1)
