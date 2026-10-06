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
res <- correlateWithBinaryPhenotype(
  traitRERw,
  fg_paths,
  min.sp         = min_sp,
  min.pos        = min_pos,
  weighted       = "auto",
  winsorizeRER   = winR,
  winsorizetrait = NULL
)
message(sprintf("[RER_BIN] Correlation done: %d genes tested.", nrow(res)))

# ── Permulation null (CC) ─────────────────────────────────────────────────────

# getPermsBinary(permmode = "cc") builds the null foreground histories with
# categoricalPermulations(): an Mk (equal-rates) transition matrix fitted on the
# observed 0/1 trait, then stochastically mapped null tip and node states, polished
# per tree. It needs no Brownian-motion simulation and no rooted tree.
#
# This requires RERconverge >= 0.3.0. Earlier builds send permmode = "cc" to
# simBinPhenoCC(), which needs an outgroup (`root_sp`) that the pipeline does not
# define, so the checks below stop the script on those builds.
num_batches     <- as.integer(args[10])
perms_per_batch <- as.integer(args[11])

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER_BIN] Permutation testing: %d batches x %d permutations (CC null) ...",
    num_batches, perms_per_batch
  ))

  if (!exists("getPermsBinary")) {
    stop("[RER_BIN] RERconverge::getPermsBinary() not found. Update RERconverge.")
  }
  if (!exists("categoricalPermulations")) {
    stop(paste0("[RER_BIN] This RERconverge build routes getPermsBinary(permmode=",
                "'cc') through simBinPhenoCC(), which needs an outgroup (root_sp) ",
                "that PhyloPhere does not define. Install RERconverge >= 0.3.0 ",
                "(provides categoricalPermulations), or set rer_perm_batches = 0."))
  }

  # getPermsBinary() scores every null with correlateWithBinaryPhenotype() at its
  # defaults: clade = "all" foreground paths, weighted = "auto", no RER winsorization,
  # min.sp = 10, min.pos = 2. permpvalcor() needs an observed reference computed the
  # same way, so res_ref uses those settings. `res` (the requested clade and
  # winsorizeRER) keeps the reported Rho and P; res_ref only serves p.perm.
  fg_paths_all <- foreground2Paths(fg_sp, geneTrees, clade = "all")
  res_ref      <- correlateWithBinaryPhenotype(traitRERw, fg_paths_all,
                                               weighted = "auto")

  run_bin_perm_batch <- function(n) getPermsBinary(
    numperms        = n,
    fg_vec          = fg_sp,
    sisters_list    = NA,          # used only when calculateenrich = TRUE
    root_sp         = NA,          # not used by the categoricalPermulations CC path
    RERmat          = traitRERw,
    trees           = geneTrees,   # the original treesObj: path length equals ncol(RERmat)
    mastertree      = geneTrees$masterTree,
    permmode        = "cc",
    method          = "k",
    calculateenrich = FALSE
  )

  message(sprintf("  [RER_BIN] Permutation batch 1 / %d", num_batches))
  perms_combined <- run_bin_perm_batch(perms_per_batch)

  if (num_batches > 1) {
    for (i in 2:num_batches) {
      message(sprintf("  [RER_BIN] Permutation batch %d / %d", i, num_batches))
      perms_combined <- combinePermData(
        perms_combined, run_bin_perm_batch(perms_per_batch), enrich = FALSE
      )
    }
  }

  n_perms <- num_batches * perms_per_batch

  # The return type of permpvalcor() depends on the RERconverge build:
  #   * bioconda v0.3.0 tag: named numeric vector with the raw proportion
  #     sum(|null| > |obs|) / N, so the (x*N + 1)/(N + 1) pseudo-count is added here.
  #   * the commit pinned in environment/install_env.sh (2bd328f7): data.frame(permpval,
  #     permstats) with a median-centered two-tailed empirical p and the
  #     (num + 1)/(denom + 1) pseudo-count already applied; permpval is used as is.
  ppc <- permpvalcor(res_ref, perms_combined)
  if (is.data.frame(ppc)) {
    permpvals <- setNames(ppc$permpval, rownames(ppc))
  } else {
    permpvals <- (as.numeric(ppc) * n_perms + 1) / (n_perms + 1)
    names(permpvals) <- names(ppc)
  }

  res$p.perm     <- permpvals[rownames(res)]
  # BH adjustment across genes
  res$p.perm.adj <- p.adjust(res$p.perm, method = "BH")
  attr(res, "n_perms") <- n_perms
  message(sprintf(
    "[RER_BIN] Permutation p-values computed for %d / %d genes (N=%d perms).",
    sum(!is.na(res$p.perm)), nrow(res), n_perms
  ))
  
  # The raw permutations object (null statistics matrices) is kept for the FCS and comparison stages
  perms_path <- sub("\\.output$", ".perms.rds", args[5])
  saveRDS(perms_combined, file = perms_path)
  message("[RER_BIN] Saved raw null permutations RDS to: ", perms_path)
} else {
  message("[RER_BIN] Permutation testing skipped (rer_perm_batches = 0).")
}

attr(res, "rer_type") <- "binary"
attr(res, "fg_sp")    <- fg_sp
attr(res, "clade")    <- rer_clade

saveRDS(res, args[5])
message("[RER_BIN] Results saved to: ", args[5])
