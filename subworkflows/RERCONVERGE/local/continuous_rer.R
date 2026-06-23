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

# ── Load trait vector and count vectors ───────────────────────────────────────
traitPath <- args[1]
load(traitPath)   # loads objects: trait_vector, n_vector, c_vector (if present)

# ── Trait transformation ───────────────────────────────────────────────────────
transform_type <- if (length(args) >= 12) args[12] else "auto"
message(sprintf("[RER] Trait transformation parameter: '%s'", transform_type))

# Check for NA values in trait_vector
vals <- trait_vector[!is.na(trait_vector)]
if (length(vals) == 0) {
  stop("ERROR: trait_vector contains no non-NA values.")
}

if (transform_type == "ha_logit") {
  if (is.null(n_vector) || is.null(c_vector)) {
    stop("ERROR: ha_logit transformation requested, but n_trait and/or c_trait were not specified or not found in the raw trait file.")
  }
  # Align species and compute Haldane-Anscombe corrected logit
  common_sp <- intersect(names(n_vector), names(c_vector))
  common_sp <- intersect(common_sp, names(trait_vector))
  if (length(common_sp) == 0) {
    stop("ERROR: No overlap between species in n_vector, c_vector, and trait_vector.")
  }
  n_val <- n_vector[common_sp]
  c_val <- c_vector[common_sp]
  
  trans_vals <- log((c_val + 0.5) / (n_val - c_val + 0.5))
  
  # Align trait_vector to these transformed values, keeping others NA
  trait_vector <- setNames(rep(NA_real_, length(trait_vector)), names(trait_vector))
  trait_vector[common_sp] <- trans_vals
  message("[RER] Applied Haldane-Anscombe corrected logit transformation")
} else if (transform_type == "logit") {
  eps  <- 1e-4
  n_clipped <- sum(vals <= 0 | vals >= 1)
  if (n_clipped > 0) {
    warning(sprintf(
      "[RER] %d value(s) outside (0,1) clipped to [%.0e, %.4f] before logit",
      n_clipped, eps, 1 - eps
    ))
  }
  trait_vector <- log(pmax(pmin(trait_vector, 1 - eps), eps) /
                      (1 - pmax(pmin(trait_vector, 1 - eps), eps)))
  message("[RER] Applied standard logit transformation")
} else if (transform_type == "arcsin") {
  # Arcsine square root: asin(sqrt(x))
  if (any(vals < 0 | vals > 1)) {
    warning("[RER] Trait values outside [0,1] detected; arcsine square root may produce NaNs.")
  }
  trait_vector <- asin(sqrt(trait_vector))
  message("[RER] Applied arcsine square root transformation")
} else if (transform_type == "log10") {
  min_nonzero  <- min(vals[vals > 0])
  shift        <- min_nonzero / 10
  trait_vector <- log10(trait_vector + shift)
  message(sprintf("[RER] Applied log10(x + %.2e) transformation", shift))
} else if (transform_type == "none") {
  message("[RER] No transformation applied (using raw values)")
} else { # auto
  lo          <- min(vals)
  hi          <- max(vals)
  unique_vals <- unique(vals)
  is_prev     <- (lo >= 0 & hi <= 1 & !all(unique_vals %in% c(0, 1)))

  if (is_prev) {
    eps          <- 1e-4
    n_clipped    <- sum(vals <= 0 | vals >= 1)
    if (n_clipped > 0) {
      warning(sprintf(
        "[RER] %d value(s) outside (0,1) clipped to [%.0e, %.4f] before logit",
        n_clipped, eps, 1 - eps
      ))
    }
    trait_vector <- log(pmax(pmin(trait_vector, 1 - eps), eps) /
                        (1 - pmax(pmin(trait_vector, 1 - eps), eps)))
    message("[RER] Auto-detected prevalence trait — applied logit transform")
  } else {
    sw <- shapiro.test(vals)
    message(sprintf(
      "[RER] Shapiro-Wilk normality test: W = %.4f, p = %.4e",
      sw$statistic, sw$p.value
    ))
    if (sw$p.value < 0.05) {
      min_nonzero  <- min(vals[vals > 0])
      shift        <- min_nonzero / 10
      trait_vector <- log10(trait_vector + shift)
      message(sprintf(
        "[RER] Auto-detected non-normal distribution — applied log10(x + %.2e) transform",
        shift
      ))
    } else {
      message("[RER] Auto-detected normal distribution — using raw values")
    }
  }
}

# ── Load gene trees ───────────────────────────────────────────────────────────
geneTrees <- readRDS(args[2])

# Resolve polytomies in the master tree: ratematrix (used by getPermsContinuous
# internally for BM simulation) requires a rooted, fully dichotomous tree.
geneTrees_di             <- geneTrees
geneTrees_di$masterTree  <- ape::multi2di(geneTrees$masterTree)

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

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER] Permutation testing: %d batches x %d permutations (BM null) ...",
    num_batches, perms_per_batch
  ))

  message(sprintf("  [RER] Permutation batch 1 / %d", num_batches))
  perms_combined <- getPermsContinuous(
    numperms        = perms_per_batch,
    traitvec        = trait_vector,
    RERmat          = traitRERw,
    annotlist       = NULL,
    trees           = geneTrees_di,
    mastertree      = geneTrees_di$masterTree,
    calculateenrich = FALSE,
    winR            = as.numeric(args[7]),
    winT            = as.numeric(args[8])
  )

  if (num_batches > 1) {
    for (i in 2:num_batches) {
      message(sprintf("  [RER] Permutation batch %d / %d", i, num_batches))
      batch_i <- getPermsContinuous(
        numperms        = perms_per_batch,
        traitvec        = trait_vector,
        RERmat          = traitRERw,
        annotlist       = NULL,
        trees           = geneTrees_di,
        mastertree      = geneTrees_di$masterTree,
        calculateenrich = FALSE,
        winR            = as.numeric(args[7]),
        winT            = as.numeric(args[8])
      )
      perms_combined <- combinePermData(perms_combined, batch_i, enrich = FALSE)
    }
  }

  permpvals <- permpvalcor(res, perms_combined)

  # Apply standard pseudo-count correction (num + 1) / (denom + 1)
  # to prevent exact 0 p-values and properly represent finite empirical probability.
  n_perms <- num_batches * perms_per_batch
  permpvals <- (permpvals * n_perms + 1) / (n_perms + 1)

  # Align by gene name (row names of res)
  res$p.perm <- permpvals[rownames(res)]
  # BH-correct the permulation p-values for multiple testing across genes.
  res$p.perm.adj <- p.adjust(res$p.perm, method = "BH")
  # Store total permutation count as attribute so the report can display
  # the minimum detectable p.perm rather than a literal zero.
  attr(res, "n_perms") <- n_perms
  message(sprintf(
    "[RER] Permutation p-values computed for %d / %d genes (N=%d perms total).",
    sum(!is.na(res$p.perm)), nrow(res), n_perms
  ))
  
  # Save the raw permutations object containing null statistics matrices for pathway-level permulations downstream
  perms_path <- sub("\\.output$", ".perms.rds", args[5])
  saveRDS(perms_combined, file = perms_path)
  message("[RER] Saved raw null permutations RDS to: ", perms_path)
} else {
  message("[RER] Permutation testing skipped (rer_perm_batches = 0).")
}

saveRDS(res, args[5])
message("[RER] Results saved to: ", args[5])
