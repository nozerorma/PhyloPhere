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

# ── Dimensional consistency guard ────────────────────────────────────────────
# foreground2Paths() derives its length from the treesObj master tree; the RER
# matrix columns were fixed against the SAME master tree at getAllResiduals()
# time. A mismatch means the foreground paths and the RER columns are recycled
# against each other during correlation — a hard error in some RERconverge
# builds, silent nonsense in others. Fail loudly.
if (length(fg_paths) != ncol(traitRERw)) {
  stop(sprintf(
    paste0("[RER_BIN] Path/RER dimension mismatch: foreground2Paths produced %d ",
           "paths but the RER matrix has %d columns. The master tree in this ",
           "treesObj is inconsistent with the one used to build the RER matrix."),
    length(fg_paths), ncol(traitRERw)
  ))
}

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
  min.sp         = min_sp,
  min.pos        = min_pos,
  weighted       = "auto",
  winsorizeRER   = winR,
  winsorizetrait = NULL
)
message(sprintf("[RER_BIN] Correlation done: %d genes tested.", nrow(res)))

# ── Permutation statistics (RERconverge binary CC permulation) ────────────────
# RERconverge::getPermsBinary(permmode = "cc") builds null foreground histories
# with categoricalPermulations(): an Mk (equal-rates) transition matrix fit on
# the observed 0/1 trait, then stochastically-mapped null tip/node states whose
# likelihood is polished per tree. No Brownian-motion simulation and no tree
# rooting are involved (that path belongs to the older simBinPhenoCC).
#
# Requires RERconverge >= 0.3.0 (categoricalPermulations()). Older builds route
# permmode = "cc" through simBinPhenoCC(), which needs a real `root_sp`
# (outgroup) that this pipeline does not define — hence the explicit version
# guard below.
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

  # getPermsBinary()'s CC branch scores every null with correlateWithBinaryPhenotype()
  # at its defaults — clade = "all" foreground paths, weighted = "auto", NO RER
  # winsorization, min.sp = 10, min.pos = 2. For permpvalcor() to compare like
  # with like, the reference observed correlation must be computed the same way.
  # `res` above (user's rer_binary_clade + winsorizeRER) still carries the
  # reported Rho / P; only the p.perm reference uses these matched settings.
  fg_paths_all <- foreground2Paths(fg_sp, geneTrees, clade = "all")
  res_ref      <- correlateWithBinaryPhenotype(traitRERw, fg_paths_all,
                                               weighted = "auto")

  run_bin_perm_batch <- function(n) getPermsBinary(
    numperms        = n,
    fg_vec          = fg_sp,
    sisters_list    = NA,          # only consumed when calculateenrich = TRUE
    root_sp         = NA,          # unused by the categoricalPermulations CC path
    RERmat          = traitRERw,
    trees           = geneTrees,   # ORIGINAL treesObj — keeps path length == ncol(RERmat)
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

  # permpvalcor()'s return type changed across RERconverge builds:
  #   * bioconda v0.3.0 tag  -> named numeric vector, raw proportion
  #     sum(|null| > |obs|) / N. Needs the (x*N + 1)/(N + 1) pseudo-count here.
  #   * install_env.sh pin (2bd328f7) -> data.frame(permpval, permstats), a
  #     median-centred two-tailed empirical p with the (num + 1)/(denom + 1)
  #     pseudo-count ALREADY applied internally. Take permpval as-is.
  ppc <- permpvalcor(res_ref, perms_combined)
  if (is.data.frame(ppc)) {
    permpvals <- setNames(ppc$permpval, rownames(ppc))
  } else {
    permpvals <- (as.numeric(ppc) * n_perms + 1) / (n_perms + 1)
    names(permpvals) <- names(ppc)
  }

  res$p.perm     <- permpvals[rownames(res)]
  # BH-correct the permulation p-values for multiple testing across genes.
  res$p.perm.adj <- p.adjust(res$p.perm, method = "BH")
  attr(res, "n_perms") <- n_perms
  message(sprintf(
    "[RER_BIN] Permutation p-values computed for %d / %d genes (N=%d perms).",
    sum(!is.na(res$p.perm)), nrow(res), n_perms
  ))
  
  # Save the raw permutations object containing null statistics matrices for pathway-level permulations downstream
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
