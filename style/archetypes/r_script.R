#!/usr/bin/env Rscript
# binary_rer.R — Compute RERconverge binary trait correlations for a given trait vector.
# PhyloPhere | subworkflows/RERCONVERGE/local/
# =============================================================================
# Called by:  RER_BIN Nextflow process (rer_bin.nf → Rscript binary_rer.R ...)
#
# Args (positional, from task.script):
#   args[1]  path to polished trait RData (loads: trait_vector, named 0/1 numeric)
#   args[2]  path to master gene-trees RDS (RERconverge format)
#   args[3]  output path for foreground-paths RDS (char2path equivalent)
#   args[4]  path to RER matrix RDS
#   args[5]  output path for binary correlation RDS
#   args[6]  min.sp            — minimum species per gene tree (integer)
#   args[7]  min.pos           — minimum foreground lineages per gene (integer)
#   args[8]  winsorizeRER      — winsorization threshold; 0 → no winsorization
#   args[9]  clade             — foreground branches: "all" | "terminal" | "ancestral"
#   args[10] rer_perm_batches  — permutation batches (0 = skip permutations)
#   args[11] rer_perms_per_batch
# =============================================================================

args <- commandArgs(TRUE)

# ── Dependencies ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr)
    library(RERconverge)
})

# ── Parse arguments ───────────────────────────────────────────────────────────

trait_path        <- args[1]
master_trees_path <- args[2]
fg_paths_out      <- args[3]
rer_matrix_path   <- args[4]
corr_out          <- args[5]
min_sp            <- as.integer(args[6])
min_pos           <- as.integer(args[7])
winsorize_rer     <- as.numeric(args[8])
clade             <- args[9]
perm_batches      <- as.integer(args[10])
perms_per_batch   <- as.integer(args[11])

# ── Load inputs ───────────────────────────────────────────────────────────────

load(trait_path)          # → trait_vector (named numeric 0/1)
master_trees <- readRDS(master_trees_path)
rer_matrix   <- readRDS(rer_matrix_path)

# ── Foreground paths ──────────────────────────────────────────────────────────

fg_paths <- char2Paths(trait_vector, master_trees, charType = clade)
saveRDS(fg_paths, file = fg_paths_out)

# ── Binary correlation ────────────────────────────────────────────────────────

winsor <- if (winsorize_rer == 0) NULL else winsorize_rer

corr <- correlateWithBinaryPhenotype(
    RERmat   = rer_matrix,
    binarytraitvec = fg_paths,
    min.sp   = min_sp,
    min.pos  = min_pos,
    winsorizeRER = winsor,
)

# ── Permutation test (optional) ───────────────────────────────────────────────

if (perm_batches > 0) {
    perm_results <- permpvalcor(
        corr,
        fg_paths,
        rer_matrix,
        trees = master_trees,
        winsorizeRER = winsor,
        numperms = perm_batches * perms_per_batch,
    )
    corr <- cbind(corr, perm_pval = perm_results$pval)
}

# ── Save output ───────────────────────────────────────────────────────────────

saveRDS(corr, file = corr_out)
message("binary_rer.R: done — ", nrow(corr), " genes")
