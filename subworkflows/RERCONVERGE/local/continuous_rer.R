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
# File: continuous_rer.R
#
# Arguments
# ---------
#   args[1]  path to trait .polished.output (RData with trait_vector)
#   args[2]  path to master gene trees RDS
#   args[3]  output path for char2path RDS
#   args[4]  path to RER matrix RDS
#   args[5]  output path for continuous correlation RDS
#   args[6]  min.sp  — minimum species per gene (integer)
#   args[7]  winsorizeRER   — winsorization threshold for RER values
#   args[8]  winsorizeTrait — winsorization threshold for trait values
#   args[9]  rer_perm_batches   — number of permutation batches (0 = skip)
#   args[10] rer_perms_per_batch — permutations per batch
#   args[11] rer_perm_mode — permutation mode: "cc" (Complete Case) or "auto"

args <- commandArgs(TRUE)

library(dplyr)
library(RERconverge)

# ── Load trait vector ─────────────────────────────────────────────────────────
traitPath <- args[1]
load(traitPath)   # loads object: trait_vector

# ── Shapiro-Wilk normality test → log10 transform if non-normal ───────────────
# As in Valenzuela et al. (2024): log10-transform traits that reject normality
# (Shapiro-Wilk p < 0.05); use raw values otherwise.
sw <- shapiro.test(trait_vector)
message(sprintf(
  "[RER] Shapiro-Wilk normality test: W = %.4f, p = %.4e",
  sw$statistic, sw$p.value
))

if (sw$p.value < 0.05) {
  min_nonzero <- min(trait_vector[trait_vector > 0], na.rm = TRUE)
  shift       <- min_nonzero / 10
  trait_vector <- log10(trait_vector + shift)
  message(sprintf(
    "[RER] Non-normal distribution — log10(x + %.2e) transform applied",
    shift
  ))
} else {
  message("[RER] Normal distribution — no transform applied, using raw values")
}

# ── Load gene trees ───────────────────────────────────────────────────────────
geneTrees <- readRDS(args[2])

# ── Convert trait vector to phylogenetic paths ────────────────────────────────
charpaths <- char2Paths(trait_vector, geneTrees)
saveRDS(charpaths, args[3])

# ── Load RER matrix ───────────────────────────────────────────────────────────
traitRERw <- readRDS(args[4])

# ── Continuous RER correlation ────────────────────────────────────────────────
message("[RER] Running correlateWithContinuousPhenotype ...")
res <- correlateWithContinuousPhenotype(
  traitRERw,
  charpaths,
  min.sp         = as.numeric(args[6]),
  winsorizeRER   = as.numeric(args[7]),
  winsorizetrait = as.numeric(args[8])
)
message(sprintf("[RER] Correlation done: %d genes tested.", nrow(res)))

# ── Permutation statistics (Brownian Motion null phenotypes) ──────────────────
# Following Valenzuela et al. (2024): 10 batches x 100 permutations per trait
# (total 1,000 permutations) using the Complete Case (CC) method.
# Empirical p-value (p.perm) = proportion of null correlations >= observed |Rho|.
num_batches      <- as.integer(args[9])
perms_per_batch  <- as.integer(args[10])
perm_mode        <- args[11]

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER] Permutation testing: %d batches x %d permutations (mode='%s') ...",
    num_batches, perms_per_batch, perm_mode
  ))

  all_perms <- vector("list", num_batches)
  for (i in seq_len(num_batches)) {
    message(sprintf("  [RER] Permutation batch %d / %d", i, num_batches))
    all_perms[[i]] <- getPermsContinuous(
      numperms   = perms_per_batch,
      traitname  = "trait",
      mastertree = geneTrees$masterTree,
      permmode   = perm_mode
    )
  }

  perms_combined <- combinePermData(all_perms)
  permpvals      <- permpvalcor(res, perms_combined)

  # Align by gene name (row names of res)
  res$p.perm <- permpvals[rownames(res)]
  message(sprintf(
    "[RER] Permutation p-values computed for %d / %d genes.",
    sum(!is.na(res$p.perm)), nrow(res)
  ))
} else {
  message("[RER] Permutation testing skipped (rer_perm_batches = 0).")
}

saveRDS(res, args[5])
message("[RER] Results saved to: ", args[5])
