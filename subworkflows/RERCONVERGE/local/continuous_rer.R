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
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "rer_permulations.R"))
source(file.path(dirname(normalizePath(script_file)), "rer_batches.R"))
source(file.path(dirname(normalizePath(script_file)), "rer_contract.R"))
batch_seed <- rer_integer(if(length(args)>=13L) args[13] else 1998, 'batch seed')
set.seed(batch_seed)
trait <- rer_trait(args[1])
transformation <- rer_transform(trait, if(length(args)>=12L) args[12] else 'auto')
trait_vector <- transformation$values
message('[RER] Effective transformation: ', transformation$effective)
min_sp <- rer_integer(args[6], 'min.sp', 3)
winR <- rer_integer(args[7], 'winsorizeRER')
winT <- rer_integer(args[8], 'winsorizeTrait')

# ── Load gene trees ───────────────────────────────────────────────────────────

geneTrees <- readRDS(args[2])

# ── Convert trait vector to phylogenetic paths ────────────────────────────────

charpaths <- char2Paths(trait_vector, geneTrees)
saveRDS(charpaths, args[3])

# ── Load RER matrix ───────────────────────────────────────────────────────────

traitRERw <- readRDS(args[4])
rer_check_matrix(traitRERw, geneTrees, names(trait_vector))

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
res <- rer_score_paths(
  traitRERw,
  charpaths,
  mode   = "continuous",
  min.sp = min_sp,
  winR   = winR,
  winT   = winT
)
res <- rer_results(res)
message(sprintf("[RER] Correlation done: %d genes tested.", nrow(res)))
attr(res, "rer_settings") <- list(mode = "continuous", min.sp = min_sp,
                                   winR = rer_winsor_limit(winR), winT = rer_winsor_limit(winT),
                                   transformation = transformation[c("requested","effective","shift")],
                                   matrix_contract = attr(traitRERw,"rer_contract"),
                                   trait_metadata = trait$metadata)

# ── Permulation null (Brownian-motion phenotypes) ─────────────────────────────

# Null phenotypes are simulated under Brownian motion (permulations, Saputra et al. 2021, who typically
# use 1,000); the configured default (conf/rerconverge.config) is 10 batches x 100 permutations.
# Gene p values use native permpvalcor on the full cycle set of each gene batch.
num_batches      <- rer_integer(args[9], "perm batches")
perms_per_batch  <- rer_integer(args[10], "perms per batch", 1)

if (num_batches > 0 && perms_per_batch > 0) {
  message(sprintf(
    "[RER] Permutation testing: %d batches x %d permutations (BM null) ...",
    num_batches, perms_per_batch
  ))

  # ── Master tree for the simulation ─────────────────────────────────────────

  # RERconverge::simpermvec uses geiger::ratematrix() and geiger::sim.char(), which
  # need a rooted, fully dichotomous tree whose tips are
  # the species of the complete-case trait vector. That tree is built here and passed
  # only as `mastertree=`. `trees=` must stay the original treesObj: replacing its
  # master tree with a multi2di() tree leaves its cached $paths, $matIndex and $ap
  # slots out of step with the topology, and char2Paths() in the null loop would
  # return paths of the wrong length (see the dimension check above).
  sim_sp     <- intersect(names(trait_vector)[!is.na(trait_vector)],
                          geneTrees$masterTree$tip.label)
  sim_trait  <- trait_vector[sim_sp]
  sim_master <- ape::keep.tip(geneTrees$masterTree, sim_sp)
  midpoint_used <- !ape::is.rooted(sim_master)
  polytomies_resolved <- !ape::is.binary(sim_master)
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

  attr(res, "rer_settings")$simulation_tree <- list(midpoint_rooted=midpoint_used,
    polytomies_resolved=polytomies_resolved, nonpositive_edges_replaced=sum(zero_edge),
    sha256=rer_hash(sim_master))
  res <- rer_write_batches(res, num_batches, perms_per_batch, batch_seed,
    function(n) rer_perms_continuous(n, sim_trait, traitRERw, geneTrees, sim_master,
                                    min.sp=min_sp, winR=winR, winT=winT),
    sub("\\.output$", ".perms.rds", args[5]))
} else {
  message("[RER] Permutation testing skipped (rer_perm_batches = 0).")
}

attr(res, "rer_type") <- "continuous"
attr(res, "rer_input_hashes") <- lapply(args[c(1,2,4)], rer_file_hash)
invisible(rer_results(res))
saveRDS(res, args[5])
message("[RER] Results saved to: ", args[5])
