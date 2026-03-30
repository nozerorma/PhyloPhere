#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: binary_rer.R
#
# Arguments
# ---------
#   args[1]  path to trait .polished.output (RData with trait_vector; 0/1 encoded)
#   args[2]  path to master gene trees RDS
#   args[3]  output path for foreground paths RDS  (char2path equivalent)
#   args[4]  path to RER matrix RDS
#   args[5]  output path for binary correlation RDS
#   args[6]  min.sp     — minimum species per gene tree (integer)
#   args[7]  min.pos    — minimum independent foreground lineages per gene (integer)
#   args[8]  winsorizeRER — winsorization threshold for RER values (0 = off → NULL)
#   args[9]  clade      — which branches to mark foreground: "all", "terminal", "ancestral"
#   args[10] rer_perm_batches   — number of permutation batches (0 = skip)
#   args[11] rer_perms_per_batch — permutations per batch

args <- commandArgs(TRUE)

library(dplyr)
library(RERconverge)

# ── Load trait vector ─────────────────────────────────────────────────────────
traitPath <- args[1]
load(traitPath)   # loads object: trait_vector  (named 0/1 numeric vector)

fg_sp <- names(trait_vector[!is.na(trait_vector) & trait_vector == 1])
bg_sp <- names(trait_vector[!is.na(trait_vector) & trait_vector == 0])
message(sprintf(
  "[RER_BIN] Trait loaded: %d foreground (1) / %d background (0) species.",
  length(fg_sp), length(bg_sp)
))
if (length(fg_sp) < 3) {
  stop("[RER_BIN] Too few foreground species (< 3). Cannot run binary analysis.")
}

# ── Load gene trees ───────────────────────────────────────────────────────────
geneTrees <- readRDS(args[2])
# Resolve polytomies in the master tree: getAllPermsBinary may need a rooted,
# fully dichotomous tree for proper permutation sampling. Dichotomize for safety.
geneTrees_di             <- geneTrees
geneTrees_di$masterTree  <- ape::multi2di(geneTrees$masterTree)
# ── Build foreground paths (binary equivalent of char2Paths) ─────────────────
# foreground2Paths assigns branch-state weights using the master tree topology.
# clade parameter controls which branches get foreground weight:
#   "all"       – transition branch + all daughter branches (broadest signal)
#   "ancestral" – only the inferred transition branch (convergence-focused)
#   "terminal"  – only terminal branches leading to foreground species
rer_clade <- args[9]
message(sprintf("[RER_BIN] Building foreground paths (clade = '%s') ...", rer_clade))
fg_paths <- foreground2Paths(fg_sp, geneTrees, clade = rer_clade)
saveRDS(fg_paths, args[3])
message(sprintf("[RER_BIN] Foreground paths computed for %d branches.", length(fg_paths)))

# ── Load RER matrix ───────────────────────────────────────────────────────────
traitRERw <- readRDS(args[4])

# ── Binary RER correlation ─────────────────────────────────────────────────────
# Uses Kendall rank correlation (unweighted when clade = "all" / "terminal" with
# 0/1 branch lengths; weighted when fractional branch lengths are present).
# winsorizeRER: pull the most extreme N RER values toward the (N+1)-th most
#   extreme before correlating — mitigates leverage from outlier branches.
#   Set to NULL (or 0 here) to skip.
winR_raw <- as.numeric(args[8])
winR     <- if (!is.na(winR_raw) && winR_raw > 0) winR_raw else NULL

min_sp  <- as.integer(args[6])
min_pos <- as.integer(args[7])

message(sprintf(
  "[RER_BIN] Running correlateWithBinaryPhenotype (min.sp=%d, min.pos=%d, winsorizeRER=%s) ...",
  min_sp, min_pos, if (is.null(winR)) "NULL" else winR
))
res <- correlateWithBinaryPhenotype(
  traitRERw,
  fg_paths,
  min.sp        = min_sp,
  min.pos       = min_pos,
  weighted      = "auto",
  winsorizeRER  = winR,
  winsorizeTrait = NULL
)
message(sprintf("[RER_BIN] Correlation done: %d genes tested.", nrow(res)))

# ── Permutation statistics (random foreground relabeling null) ────────────────
# Binary permutations randomly reassign foreground labels among species while
# preserving the count of foreground species, providing an empirical null.
# Note: unlike continuous (BM simulation), this uses permutations of labels.
num_batches     <- as.integer(args[10])
perms_per_batch <- as.integer(args[11])

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER_BIN] Permutation testing: %d batches x %d permutations ...",
    num_batches, perms_per_batch
  ))

  message(sprintf("  [RER_BIN] Permutation batch 1 / %d", num_batches))
  perms_combined <- getAllPermsBinary(
    numperms   = perms_per_batch,
    fg_sp      = fg_sp,
    RERmat     = traitRERw,
    trees      = geneTrees_di,
    mastertree = geneTrees_di$masterTree
  )

  if (num_batches > 1) {
    for (i in 2:num_batches) {
      message(sprintf("  [RER_BIN] Permutation batch %d / %d", i, num_batches))
      batch_i <- getAllPermsBinary(
        numperms   = perms_per_batch,
        fg_sp      = fg_sp,
        RERmat     = traitRERw,
        trees      = geneTrees_di,
        mastertree = geneTrees_di$masterTree
      )
      perms_combined <- combinePermData(perms_combined, batch_i, enrich = FALSE)
    }
  }

  permpvals      <- permpvalcor(res, perms_combined)
  res$p.perm     <- permpvals[rownames(res)]
  attr(res, "n_perms") <- num_batches * perms_per_batch
  message(sprintf(
    "[RER_BIN] Permutation p-values computed for %d / %d genes (N=%d perms).",
    sum(!is.na(res$p.perm)), nrow(res), num_batches * perms_per_batch
  ))
} else {
  message("[RER_BIN] Permutation testing skipped (rer_perm_batches = 0).")
}

attr(res, "rer_type") <- "binary"
attr(res, "fg_sp")    <- fg_sp
attr(res, "clade")    <- rer_clade

saveRDS(res, args[5])
message("[RER_BIN] Results saved to: ", args[5])
