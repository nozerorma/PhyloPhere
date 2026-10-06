# continuous_rer.R — RERconverge correlation of gene RERs with a continuous trait.
# PhyloPhere | subworkflows/RERCONVERGE/local/
# =============================================================================
# Called by:  RER_CONT Nextflow process (rer_cont.nf → Rscript continuous_rer.R ...)
#
# Transforms the trait (args[12]), converts it to phylogenetic paths on the master tree,
# correlates the paths with the RER matrix (correlateWithContinuousPhenotype) and,
# when permutations are requested, adds p.perm (empirical, from Brownian-motion null
# phenotypes) and p.perm.adj (BH). The result carries the attribute n_perms.
#
# Args (positional, from task.script):
#   args[1]  polished trait RData (trait_vector, n_vector, c_vector)
#   args[2]  master gene-trees RDS (RER_TREES)
#   args[3]  output: char2Paths RDS
#   args[4]  RER matrix RDS (RER_MATRIX)
#   args[5]  output: continuous correlation RDS (the permutations go to the same name
#            with .output replaced by .perms.rds)
#   args[6]  min.sp: minimum species per gene (integer)
#   args[7]  winsorizeRER: winsorization threshold of the RER values
#   args[8]  winsorizeTrait: winsorization threshold of the trait values
#   args[9]  rer_perm_batches: permutation batches (0 = no permutations)
#   args[10] rer_perms_per_batch: permutations per batch
#   args[11] permutation mode ("cc"); passed by the process, not read by this script
#   args[12] trait transformation: auto | ha_logit | logit | arcsin | log10 | none
# =============================================================================

args <- commandArgs(TRUE)

# ── Dependencies ──────────────────────────────────────────────────────────────

library(dplyr)
library(RERconverge)

# ── Load trait vector and count vectors ───────────────────────────────────────

traitPath <- args[1]
load(traitPath)   # trait_vector, n_vector, c_vector (n_vector and c_vector may be NULL)

# ── Trait transformation ───────────────────────────────────────────────────────

transform_type <- if (length(args) >= 12) args[12] else "auto"
message(sprintf("[RER] Trait transformation parameter: '%s'", transform_type))

# At least one non-NA value is needed
vals <- trait_vector[!is.na(trait_vector)]
if (length(vals) == 0) {
  stop("ERROR: trait_vector contains no non-NA values.")
}

if (transform_type == "ha_logit") {
  if (is.null(n_vector) || is.null(c_vector)) {
    stop("ERROR: ha_logit transformation requested, but n_trait and/or c_trait were not specified or not found in the raw trait file.")
  }
  # Haldane-Anscombe logit of c_trait out of n_trait, on the species present in all three vectors
  common_sp <- intersect(names(n_vector), names(c_vector))
  common_sp <- intersect(common_sp, names(trait_vector))
  if (length(common_sp) == 0) {
    stop("ERROR: No overlap between species in n_vector, c_vector, and trait_vector.")
  }
  n_val <- n_vector[common_sp]
  c_val <- c_vector[common_sp]
  
  trans_vals <- log((c_val + 0.5) / (n_val - c_val + 0.5))
  
  # The transformed values replace trait_vector; the other species become NA
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
  # asin(sqrt(x)), defined for proportions in [0, 1]
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

# ── Convert trait vector to phylogenetic paths ────────────────────────────────

charpaths <- char2Paths(trait_vector, geneTrees)
saveRDS(charpaths, args[3])

# ── Load RER matrix ───────────────────────────────────────────────────────────

traitRERw <- readRDS(args[4])

# ── Dimension check ───────────────────────────────────────────────────────────

# The number of paths from char2Paths() comes from allPaths(treesObj$masterTree), and the
# RER matrix columns come from the master tree used in getAllResiduals(). If they
# differ, the trait paths are recycled against the RER columns: an error inside
# getAllCor() ("logical subscript too long") or, depending on the RERconverge build,
# meaningless correlations. The script stops instead.
if (length(charpaths) != ncol(traitRERw)) {
  stop(sprintf(
    paste0("[RER] Path/RER dimension mismatch: char2Paths produced %d paths but ",
           "the RER matrix has %d columns. The master tree in this treesObj is ",
           "inconsistent with the one used to build the RER matrix."),
    length(charpaths), ncol(traitRERw)
  ))
}

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

# ── Permulation null (Brownian-motion phenotypes) ─────────────────────────────

# Null phenotypes are simulated under Brownian motion; the configured default (conf/
# rerconverge.config) is 10 batches x 100 permutations, as in Valenzuela et al. (2024).
# p.perm is the empirical p-value computed by permpvalcor() from the null correlations.
num_batches      <- as.integer(args[9])
perms_per_batch  <- as.integer(args[10])

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER] Permutation testing: %d batches x %d permutations (BM null) ...",
    num_batches, perms_per_batch
  ))

  # ── Master tree for the simulation ─────────────────────────────────────────

  # getPermsContinuous() simulates the null phenotypes with geiger::ratematrix() and
  # geiger::sim.char(), which need a rooted, fully dichotomous tree whose tips are
  # the species of the complete-case trait vector. That tree is built here and passed
  # only as `mastertree=`. `trees=` must stay the original treesObj: replacing its
  # master tree with a multi2di() tree leaves its cached $paths, $matIndex and $ap
  # slots out of step with the topology, and char2Paths() in the null loop would
  # return paths of the wrong length (see the dimension check above).
  sim_sp     <- intersect(names(trait_vector)[!is.na(trait_vector)],
                          geneTrees$masterTree$tip.label)
  sim_trait  <- trait_vector[sim_sp]
  sim_master <- ape::keep.tip(geneTrees$masterTree, sim_sp)
  if (!ape::is.rooted(sim_master)) {
    # The pipeline defines no outgroup, so the tree is midpoint-rooted; the Brownian-motion
    # simulation is only weakly sensitive to the placement of the root.
    sim_master <- phytools::midpoint.root(sim_master)
  }
  sim_master <- ape::multi2di(sim_master)
  # multi2di() adds zero-length edges; they are set to 1e-8 to keep ratematrix()
  # non-singular.
  zero_edge <- sim_master$edge.length <= 0
  if (any(zero_edge)) {
    sim_master$edge.length[zero_edge] <- 1e-8
  }
  message(sprintf(
    "[RER] BM-null master tree: %d tips (rooted=%s, binary=%s)",
    length(sim_master$tip.label),
    ape::is.rooted(sim_master), ape::is.binary(sim_master)
  ))

  message(sprintf("  [RER] Permutation batch 1 / %d", num_batches))
  perms_combined <- getPermsContinuous(
    numperms        = perms_per_batch,
    traitvec        = sim_trait,
    RERmat          = traitRERw,
    annotlist       = NULL,
    trees           = geneTrees,
    mastertree      = sim_master,
    calculateenrich = FALSE,
    winR            = as.numeric(args[7]),
    winT            = as.numeric(args[8])
  )

  if (num_batches > 1) {
    for (i in 2:num_batches) {
      message(sprintf("  [RER] Permutation batch %d / %d", i, num_batches))
      batch_i <- getPermsContinuous(
        numperms        = perms_per_batch,
        traitvec        = sim_trait,
        RERmat          = traitRERw,
        annotlist       = NULL,
        trees           = geneTrees,
        mastertree      = sim_master,
        calculateenrich = FALSE,
        winR            = as.numeric(args[7]),
        winT            = as.numeric(args[8])
      )
      perms_combined <- combinePermData(perms_combined, batch_i, enrich = FALSE)
    }
  }

  n_perms <- num_batches * perms_per_batch

  # The return type of permpvalcor() depends on the RERconverge build:
  #   * bioconda v0.3.0 tag: named numeric vector with the raw proportion
  #     sum(|null| > |obs|) / N, so the (x*N + 1)/(N + 1) pseudo-count is added here.
  #   * the commit pinned in environment/install_env.sh (2bd328f7): data.frame(permpval,
  #     permstats) with a median-centered two-tailed empirical p and the
  #     (num + 1)/(denom + 1) pseudo-count already applied; permpval is used as is.
  ppc <- permpvalcor(res, perms_combined)
  if (is.data.frame(ppc)) {
    permpvals <- setNames(ppc$permpval, rownames(ppc))
  } else {
    permpvals <- (as.numeric(ppc) * n_perms + 1) / (n_perms + 1)
    names(permpvals) <- names(ppc)
  }

  # Match by gene name (row names of res)
  res$p.perm <- permpvals[rownames(res)]
  # BH adjustment across genes
  res$p.perm.adj <- p.adjust(res$p.perm, method = "BH")
  # The total number of permutations lets the report show the smallest detectable
  # p.perm instead of 0.
  attr(res, "n_perms") <- n_perms
  message(sprintf(
    "[RER] Permutation p-values computed for %d / %d genes (N=%d perms total).",
    sum(!is.na(res$p.perm)), nrow(res), n_perms
  ))
  
  # The raw permutations object (null statistics matrices) is kept for the FCS and comparison stages
  perms_path <- sub("\\.output$", ".perms.rds", args[5])
  saveRDS(perms_combined, file = perms_path)
  message("[RER] Saved raw null permutations RDS to: ", perms_path)
} else {
  message("[RER] Permutation testing skipped (rer_perm_batches = 0).")
}

saveRDS(res, args[5])
message("[RER] Results saved to: ", args[5])
