# binary_rer.R — RERconverge correlation of gene RERs with a binary (0/1) trait.
# PhyloPhere | subworkflows/RERCONVERGE/local/
# =============================================================================
# Called by:  RER_BIN Nextflow process (rer_bin.nf → Rscript binary_rer.R ...)
#
# Builds the foreground paths of the 0/1 trait on the master tree, correlates them with
# the RER matrix (correlateWithBinaryPhenotype) and, when permutations are requested,
# adds an empirical p-value (p.perm) and its BH adjustment (p.perm.adj) from the
# RERconverge categorical (CC) permulation. The result carries the attributes
# rer_type, fg_sp, clade and n_perms.
#
# Args (positional, from task.script):
#   args[1]  polished trait RData (trait_vector, named 0/1 numeric vector)
#   args[2]  master gene-trees RDS (RER_TREES)
#   args[3]  output: foreground paths RDS
#   args[4]  RER matrix RDS (RER_MATRIX)
#   args[5]  output: binary correlation RDS (the permutations go to the same name
#            with .output replaced by .perms.rds)
#   args[6]  min.sp: minimum species per gene (integer)
#   args[7]  min.pos: minimum independent foreground lineages per gene (integer)
#   args[8]  winsorizeRER: winsorization threshold of the RER values (0 = none)
#   args[9]  clade: foreground branches, "all" | "terminal" | "ancestral"
#   args[10] rer_perm_batches: permutation batches (0 = no permutations)
#   args[11] rer_perms_per_batch: permutations per batch
# =============================================================================

args <- commandArgs(TRUE)

# ── Dependencies ──────────────────────────────────────────────────────────────

library(dplyr)
library(RERconverge)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "rer_permulations.R"))
source(file.path(dirname(normalizePath(script_file)), "rer_batches.R"))
batch_seed <- if (length(args) >= 12) as.integer(args[12]) else 1998L
set.seed(batch_seed)


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

# ── Foreground paths ──────────────────────────────────────────────────────────

# foreground2Paths() weights the branches of the master tree as foreground; `clade`
# selects which ones:
#   "all"       the transition branch and all its daughter branches (broadest)
#   "ancestral" only the inferred transition branch
#   "terminal"  only the terminal branches of the foreground species
rer_clade <- args[9]
message(sprintf("[RER_BIN] Building foreground paths (clade = '%s') ...", rer_clade))
fg_paths <- foreground2Paths(fg_sp, geneTrees, clade = rer_clade)
saveRDS(fg_paths, args[3])
message(sprintf("[RER_BIN] Foreground paths computed for %d branches.", length(fg_paths)))

# ── Load RER matrix ───────────────────────────────────────────────────────────
traitRERw <- readRDS(args[4])

# ── Dimension check ───────────────────────────────────────────────────────────

# The number of foreground paths comes from the master tree of the treesObj, and the
# RER matrix columns come from the master tree used in getAllResiduals(). If they
# differ, the correlation recycles one against the other: an error in some RERconverge
# builds and meaningless values in others. The script stops instead.
if (length(fg_paths) != ncol(traitRERw)) {
  stop(sprintf(
    paste0("[RER_BIN] Path/RER dimension mismatch: foreground2Paths produced %d ",
           "paths but the RER matrix has %d columns. The master tree in this ",
           "treesObj is inconsistent with the one used to build the RER matrix."),
    length(fg_paths), ncol(traitRERw)
  ))
}

# ── Binary RER correlation ─────────────────────────────────────────────────────

# weighted = "auto" lets RERconverge choose between the unweighted test (0/1 foreground
# weights) and the weighted one (fractional weights). winsorizeRER pulls the most
# extreme RER values toward the next most extreme one before correlating, which limits
# the leverage of outlier branches; 0 or NA gives NULL (no winsorization).
winR_raw <- as.numeric(args[8])
winR     <- if (!is.na(winR_raw) && winR_raw > 0) winR_raw else NULL

min_sp  <- as.integer(args[6])
min_pos <- as.integer(args[7])

message(sprintf(
  "[RER_BIN] Running correlateWithBinaryPhenotype (min.sp=%d, min.pos=%d, winsorizeRER=%s) ...",
  min_sp, min_pos, if (is.null(winR)) "NULL" else winR
))
res <- rer_score_paths(
  traitRERw,
  fg_paths,
  mode    = "binary",
  min.sp  = min_sp,
  min.pos = min_pos,
  winR    = winR
)
message(sprintf("[RER_BIN] Correlation done: %d genes tested.", nrow(res)))
attr(res, "rer_settings") <- list(mode = "binary", min.sp = min_sp, min.pos = min_pos,
                                   winR = winR, clade = rer_clade)

# ── Permulation null (CC) ─────────────────────────────────────────────────────

# The CC null builds foreground histories with categoricalPermulations(): an Mk
# (equal-rates) transition matrix fitted on the
# observed 0/1 trait, then stochastically mapped null tip and node states, polished
# per tree. It needs no Brownian-motion simulation and no rooted tree.
#
# rer_perms_binary scores these histories with the observed correlation settings.
# categoricalPermulations is required; the legacy CC simulator needs an outgroup.
num_batches     <- as.integer(args[10])
perms_per_batch <- as.integer(args[11])

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER_BIN] Permutation testing: %d batches x %d permutations (CC null) ...",
    num_batches, perms_per_batch
  ))

  if (!exists("categoricalPermulations")) {
    stop(paste0("[RER_BIN] categoricalPermulations() not found. Install the pinned ",
                "RERconverge build providing categoricalPermulations, ",
                "or set rer_perm_batches = 0."))
  }

  run_bin_perm_batch <- function(n) rer_perms_binary(
    numperms        = n,
    fg_vec          = fg_sp,
    RERmat          = traitRERw,
    trees           = geneTrees,   # the original treesObj: path length equals ncol(RERmat)
    min.sp          = min_sp,
    min.pos         = min_pos,
    winR            = winR,
    clade           = rer_clade
  )

  res <- rer_write_batches(res, num_batches, perms_per_batch, batch_seed,
                           run_bin_perm_batch, sub("\\.output$", ".perms.rds", args[5]))
} else {
  message("[RER_BIN] Permutation testing skipped (rer_perm_batches = 0).")
}

attr(res, "rer_type") <- "binary"
attr(res, "fg_sp")    <- fg_sp
attr(res, "clade")    <- rer_clade

saveRDS(res, args[5])
message("[RER_BIN] Results saved to: ", args[5])
