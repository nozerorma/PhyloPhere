#!/usr/bin/env Rscript
# permulations.R — Harvest a pool of permulated FG/BG labelings that match the observed contrast design.
# PhyloPhere | subworkflows/CT/local/scripts/
# =============================================================================
# Called by:  RESAMPLE Nextflow process (ct_resample.nf → Rscript permulations.R ...)
#
# Args (positional, from task.script; the first six are required):
#   args[1]  tree                 species tree (Newick)
#   args[2]  config               observed traitfile (V1 species, V2 label 0/1, V3 pair id): sets
#                                 N_pairs_obs and, with multi_hypothesis, locates
#                                 contrast_hypotheses_pairs.tsv next to it
#   args[3]  cycles               size of the pool to harvest (params.caas_full_perms)
#   args[4]  strategy             auto | OU | BM (params.perm_strategy)
#   args[5]  phenotypes           trait table (species + value, optionally n and c count columns)
#   args[6]  outdir               output directory
#   args[7]  chunk_size           permulations per resample_NNN.tab (default 500)
#   args[8]  pss_top_pct          fraction of candidate pairs kept by PSS (default 0.01)
#   args[9]  max_tries            initial draw budget (default 1e6)
#   args[10] pheno_col            value column of the trait table ("" = first numeric)
#   args[11] n_col, args[12] c_col   count columns for Jeffreys intervals ("" = auto-detect)
#   args[13] resample_use_n       use the count columns when present
#   args[14] trait_type           auto | continuous | ordinal
#   args[15] multi_hypothesis     also harvest FOP hypotheses for every accepted cycle
#   args[16] max_fop              maximum FOP hypotheses (H1..Hn) per cycle (default 100)
#   args[17] n_cpus               workers of the FOP harvest (default: SLURM allocation, else cores)
#   args[18] seed                 RNG seed; required with multi_hypothesis
#   args[19] match_pss            match the PSS profile of the observed canonical pairs (default true; not applied to count traits)
#   args[20] match_pss_tol        relative PSS tolerance of that matching (default 0.25)
#
# Method: the trait is simulated under BM on the species tree (rescaled to the fitted OU
# when OU is the selected model), and the simulated ranks are mapped back onto the observed
# values, so each permulation has the observed marginal distribution. Strategy auto
# selects BM or OU by AIC, as the observed selector does.
#
# Pool harvesting: every accepted permulation carries exactly N_pairs_obs pairs (the contrast
# selection is lean_contrast_selector.R). By default (match_pss) the pairs are chosen one by one so
# that their PSS follows the PSS of the observed canonical pairs, and a permulation is accepted only
# if every pair is within match_pss_tol of its observed counterpart and all pairs are independent
# (Tier 1). Without matching, or for count traits, the pairs are the N_pairs_obs best ones, graded
# by the modified Dunn index:
#   Tier 1 : all pairs have mod_dunn >= 1 (fully independent)
#   Tier 2 : exactly one pair falls below mod_dunn 1 (only used to fill a shortfall)
#
# Outputs (outdir): resample_NNN.tab (cycle, fg, bg; no header), permulation_manifest.tsv
# (one row per cycle: tier, pair count, Dunn, trait values, design of the canonical pairs) and
# permulation_harvest.tsv (draws, rejections by reason, acceptance by PSS tolerance, pool and capacity of a
# sample of the draws), and with multi_hypothesis fop_labelings.tab and fop_pairs.tsv.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(parallel)
})

# Parallel-safe RNG for the forked workers of the FOP mirror. The generator kind is set
# unconditionally, so that it is the same with or without a seed; the seed itself is
# applied after the arguments are parsed.
RNGkind("L'Ecuyer-CMRG")

log_msg <- function(tag, ...) write(paste0("[", tag, "] ", format(Sys.time()), " ", paste0(...)), stdout())

# The evolutionary model fit, the BM/OU covariances and the AIC model selection come from
# the vendored phyloq engine (pss_core.R), sourced below once script_dir is known.

# ── Permulation primitives (RERconverge) ──────────────────────────────────────

# One BM simulation of the trait on the tree, with the rate matrix estimated from `namedvec`.
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm   <- ratematrix(treewithbranchlengths, namedvec)
  sims <- sim.char(treewithbranchlengths, rm, nsim = 1)
  setNames(as.data.frame(sims)[, 1], rownames(sims))
}

# Rank matching: the simulated vector supplies the ordering, the observed vector supplies
# the values. The permulated marginal distribution is therefore identical to the observed
# one by construction.
simpermvec <- function(namedvec, treewithbranchlengths) {
  vec       <- simulatevec(namedvec, treewithbranchlengths)
  simsorted <- sort(vec)
  simsorted[] <- sort(namedvec)
  simsorted
}

# ── CLI ───────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: permulations.R <tree> <config> <cycles> <strategy> <phenotypes> <outdir> ",
       "[chunk_size] [pss_top_pct] [max_tries] [pheno_col] ",
       "[n_col] [c_col] [resample_use_n] [trait_type] [multi_hypothesis] [max_fop] ",
       "[n_cpus] [seed]")
}

arg_or <- function(i, default, cast = as.character) {
  if (length(args) >= i && nzchar(args[i])) cast(args[i]) else default
}

tree.path          <- args[1]
config.file        <- args[2]
number.of.cycles   <- as.integer(args[3])   # size of the harvested pool
selection.strategy <- tolower(args[4])      # "ou" | "bm"
phenotypes         <- args[5]
outdir             <- args[6]
chunk.size         <- arg_or(7,  500L,       as.integer)
pss_top_pct        <- arg_or(8,  0.01,       as.numeric)
max_tries          <- arg_or(9,  1000000L,   as.integer)
pheno_col_name     <- arg_or(10, "")
n_col              <- arg_or(11, "")   # denominator column (e.g. adult_necropsy_count)
c_col              <- arg_or(12, "")   # numerator column   (e.g. malignant_count)
resample_use_n     <- tolower(arg_or(13, "true")) %in% c("1", "true", "t", "yes", "y")
trait_type         <- tolower(arg_or(14, "auto"))
# The null mirrors the observed design: with multi_hypothesis every accepted cycle also gets a FOP
# hypothesis harvest (fop_labelings.tab, fop_pairs.tsv), so that it is pooled like the observed data.
fop_null           <- tolower(arg_or(15, "false")) %in% c("1", "true", "t", "yes", "y")
max_fop            <- arg_or(16, 100L, as.integer)

# ── Parallelism and RNG seed ──────────────────────────────────────────────────
# n_cpus drives the forked FOP-mirror harvest only (the pool harvest is serial). Without
# the argument it falls back to the SLURM allocation, then to the detected cores, and it
# never exceeds the cores that are available.
.detected_cores <- tryCatch(parallel::detectCores(), error = function(e) 1L)
if (!is.finite(.detected_cores) || .detected_cores < 1L) .detected_cores <- 1L
n_cpus <- arg_or(17, NA_integer_, as.integer)
if (is.na(n_cpus)) {
  .slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  n_cpus <- if (!is.na(.slurm_cpus) && .slurm_cpus >= 1L) .slurm_cpus else .detected_cores
}
n_cpus <- max(1L, min(as.integer(n_cpus), .detected_cores))

# With a seed (the pipeline passes params.seed, 1998 by default) the parallel streams are
# reproducible; without one the RNG is left unseeded.
seed_arg <- arg_or(18, NA_integer_, as.integer)
match_pss     <- tolower(arg_or(19, "true")) %in% c("1", "true", "t", "yes", "y")
match_pss_tol <- arg_or(20, 0.25, as.numeric)
if (!is.finite(match_pss_tol) || match_pss_tol <= 0) stop("match_pss_tol must be a positive number.")
if (!is.na(seed_arg)) {
  set.seed(seed_arg)
  log_msg("INFO", sprintf("RNG seeded with %d (L'Ecuyer-CMRG); FOP mirror uses %d core(s)",
                          seed_arg, n_cpus))
} else {
  log_msg("INFO", sprintf("RNG unseeded; FOP mirror uses %d core(s)", n_cpus))
}

if (fop_null && is.na(seed_arg)) {
  stop("permulations.R: the FOP harvest needs the pipeline seed (argument 18, params.seed).")
}

# Design matching keeps only null cycles whose FOP harvest yields at least as many
# hypotheses as the observed harvest (see "Design matching" below). It is requested here
# and takes effect only when the FOP harvest is on and the observed count is known.
match_fop <- TRUE

if (!selection.strategy %in% c("auto", "best_model", "ou", "bm")) {
  log_msg("WARN", sprintf("Unknown strategy '%s', defaulting to 'auto'", selection.strategy))
  selection.strategy <- "auto"
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  log_msg("INFO", "Created output directory: ", outdir)
}

# ── Locate the selector and the PSS engine next to this script ────────────────
script_dir <- {
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]))) else getwd()
}
lean_script <- file.path(script_dir, "lean_contrast_selector.R")

# ── Config: species | binary label | pair (cluster) id ────────────────────────
cfg <- read.table(config.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
foreground.species <- cfg$V1[cfg$V2 == "1"]
background.species <- cfg$V1[cfg$V2 == "0"]

if (ncol(cfg) < 3) {
  stop("Discovery config '", config.file, "' has no V3 cluster-id column, so the ",
       "observed independent pair count cannot be determined. Permulations must ",
       "match the pair count produced by the contrast selection step.")
}
target_pairs <- suppressWarnings(max(as.integer(cfg$V3), na.rm = TRUE))
if (!is.finite(target_pairs) || target_pairs <= 0L) {
  stop("Could not parse a positive pair count from the config V3 column.")
}
n_fg <- sum(cfg$V2 == "1"); n_bg <- sum(cfg$V2 == "0")
if (n_fg != target_pairs || n_bg != target_pairs) {
  log_msg("WARN", sprintf(
    "Config has %d FG / %d BG species but %d pair ids — permulations will target %d pairs",
    n_fg, n_bg, target_pairs, target_pairs))
}
log_msg("INFO", sprintf("Observed independent pair count from config V3: N_pairs_obs = %d", target_pairs))

# Observed FOP hypothesis count: the number of hypotheses the observed harvest produced
# (at most max_fop), read from contrast_hypotheses_pairs.tsv next to the discovery config.
# It is the design size to which every null cycle is matched.
n_hyp_obs <- 0L
if (fop_null) {
  .cfg_dir <- if (dir.exists(config.file)) config.file else dirname(config.file)
  .hp_file <- file.path(.cfg_dir, "contrast_hypotheses_pairs.tsv")
  if (file.exists(.hp_file)) {
    .hp <- read.delim(.hp_file, stringsAsFactors = FALSE)
    if ("hypothesis_id" %in% names(.hp)) n_hyp_obs <- length(unique(.hp$hypothesis_id))
  }
  if (n_hyp_obs > 0L) {
    log_msg("INFO", sprintf("Observed FOP harvest: %d hypotheses; null cycles must reach >= %d", n_hyp_obs, n_hyp_obs))
  } else {
    log_msg("WARN", "no observed contrast_hypotheses_pairs.tsv next to the config; null cycles are not design-matched")
  }
}
match_fop <- fop_null && match_fop && n_hyp_obs > 0L

first_line <- readLines(phenotypes, n = 1L, warn = FALSE)
delim_char <- if (grepl(",", first_line) && !grepl("\t", first_line)) "," else "\t"
first_flds <- strsplit(first_line, delim_char, fixed = TRUE)[[1]]
has_header <- length(first_flds) < 2 || is.na(suppressWarnings(as.numeric(first_flds[2])))

pheno_raw <- read.delim(phenotypes, sep = delim_char, header = has_header,
                        stringsAsFactors = FALSE,
                        na.strings = c("", "NA", "NaN", "nan", "NULL", "null"))
log_msg("INFO", sprintf("Read %d phenotype rows from %s (header=%s)",
                        nrow(pheno_raw), basename(phenotypes), has_header))

sp_col <- if ("species" %in% names(pheno_raw)) "species" else names(pheno_raw)[1]
val_col <- if (nzchar(pheno_col_name) && pheno_col_name %in% names(pheno_raw)) {
  pheno_col_name
} else {
  cand <- setdiff(names(pheno_raw)[vapply(pheno_raw, is.numeric, logical(1))], sp_col)
  if (length(cand)) cand[1] else names(pheno_raw)[2]
}
log_msg("INFO", sprintf("Using species column '%s' and value column '%s'", sp_col, val_col))

phenotype.df <- data.frame(
  species = as.character(pheno_raw[[sp_col]]),
  value   = suppressWarnings(as.numeric(pheno_raw[[val_col]])),
  stringsAsFactors = FALSE
)

# ── Count data → Jeffreys CIs ─────────────────────────────────────────────────
if (resample_use_n) {
  if (!nzchar(n_col) || !n_col %in% names(pheno_raw)) {
    cand_n <- c("n_trait", "n", "n_population", "sample_size", "total_count", "N")
    found_n <- intersect(cand_n, names(pheno_raw))
    if (length(found_n) > 0) n_col <- found_n[1]
  }
  if (!nzchar(c_col) || !c_col %in% names(pheno_raw)) {
    cand_c <- c("c_trait", "c", "c_cases", "n_cases", "cases", "C")
    found_c <- intersect(cand_c, names(pheno_raw))
    if (length(found_c) > 0) c_col <- found_c[1]
  }
}

use_ci <- resample_use_n && nzchar(n_col) && nzchar(c_col) &&
          n_col %in% names(pheno_raw) && c_col %in% names(pheno_raw)
if (use_ci) {
  if (!requireNamespace("binom", quietly = TRUE)) {
    stop("Count columns were supplied but the 'binom' package is unavailable.")
  }
  np <- suppressWarnings(as.numeric(pheno_raw[[n_col]]))
  nc <- suppressWarnings(as.numeric(pheno_raw[[c_col]]))
  good <- is.finite(np) & is.finite(nc) & np > 0 & nc >= 0 & nc <= np
  phenotype.df$n_pop <- np   # denominator, for the pair_n ranking tiebreak
  phenotype.df$ci_lb <- NA_real_; phenotype.df$ci_ub <- NA_real_
  if (any(good)) {
    b <- binom::binom.confint(nc[good], np[good], method = "bayes",
                              priors = c(0.5, 0.5), conf.level = 0.95)
    phenotype.df$ci_lb[good] <- b[, 5]; phenotype.df$ci_ub[good] <- b[, 6]
  }
  log_msg("INFO", sprintf("Count columns '%s'/'%s' found: %d/%d species have Jeffreys CIs",
                          n_col, c_col, sum(good), length(good)))
}

keep <- !is.na(phenotype.df$species) & is.finite(phenotype.df$value)
if (use_ci) keep <- keep & is.finite(phenotype.df$ci_lb) & is.finite(phenotype.df$ci_ub)
phenotype.df <- phenotype.df[keep, , drop = FALSE]
if (nrow(phenotype.df) == 0) stop("No valid phenotype rows remain after NA filtering")

# ── Vendored phyloq engine and shared selector ────────────────────────────────
pss_core_script <- file.path(script_dir, "pss_core.R")
for (.f in c(lean_script, pss_core_script)) {
  if (!file.exists(.f)) stop(basename(.f), " not found next to permulations.R (looked in '", script_dir, "').")
}
source(pss_core_script)   # fit_models / covariances_from_fits / select_model / calculate_pairwise_scores
source(lean_script)       # selection_context + shared selection core
log_msg("INFO", "Shared contrast selector + phyloq PSS loaded from ", script_dir)

# ── Selection context ─────────────────────────────────────────────────────────
# Tree, distances and evolutionary model, fitted once on the observed trait by the same
# selection_context() that the observed selector uses. The covariances stay fixed across
# all permulation draws.
.force_model <- if (selection.strategy %in% c("bm", "ou")) toupper(selection.strategy) else NULL
ctx <- selection_context(setNames(phenotype.df$value, phenotype.df$species),
                         read.tree(tree.path), force_model = .force_model)
pruned.tree     <- ctx$tree
starting.values <- ctx$trait_vec
D               <- ctx$D
obs_fits        <- ctx$fits
selected_model  <- ctx$selected_model
cov_bm          <- ctx$cov_bm
cov_ou          <- ctx$cov_ou
log_msg("INFO", sprintf("Evolutionary model: %s (AIC BM = %.3f, OU = %.3f; delta = %.3f)",
                        selected_model, fit_aic(obs_fits$BM), fit_aic(obs_fits$OU),
                        fit_aic(obs_fits$BM) - fit_aic(obs_fits$OU)))

# ── Ultrametric check (warn-only) ─────────────────────────────────────────────
# Same check as in subworkflows/TRAIT_ANALYSIS/local/src/commons.R; keep the two in sync.
# Contrast independence (Dunn) and the OU/BM PSS assume a time tree: on an ML phylogram a
# long terminal branch (rate, not time) inflates the diameter of a species' contrast pair
# and wrongly fails otherwise-independent pairs. PhyloPhere expects a dated species tree,
# so the script warns instead of rate-smoothing an arbitrary phylogram.
if (!is.ultrametric(pruned.tree, tol = 1e-6)) {
  log_msg("WARN", "permulation tree is NOT ultrametric (phylogram); contrast ",
          "independence assumes a time tree. Supply a dated species tree.")
}

# Simulation tree for the rank-match null: OU-rescaled when OU is selected.
if (selected_model == "OU") {
  simulation_tree <- rescale(pruned.tree, "OU", as.numeric(obs_fits$OU$opt$alpha))
  simulation_tree$edge.length[simulation_tree$edge.length <= 0] <- 1e-8
} else {
  simulation_tree <- pruned.tree
}

# ── PSS profile of the observed canonical pairs ───────────────────────────────
# The observed selector chooses its pairs until the Dunn index stops it, so they are the pairs the
# trait offers; the permulations are matched to their PSS (see match_pss_select). The profile is
# computed with the same candidate function on the observed trait, in the order of the pair ids of
# the observed traitfile (the order of selection). Count traits (Jeffreys CI gate) are not matched.
pss_profile <- NULL
if (match_pss && use_ci) {
  log_msg("INFO", "PSS matching is not applied to count traits (CI gate): the null keeps the best Dunn-independent pairs")
} else if (match_pss) {
  obs_cand <- lean_candidate_df(
    starting.values, D, target_pairs, pruned.tree, cov_bm, cov_ou, selected_model,
    top_pct = 1,
    ordinal = if (trait_type == "ordinal") TRUE else if (trait_type == "continuous") FALSE else NULL)$cand_df
  pid <- sort(unique(as.integer(cfg$V3)))
  obs_pss <- vapply(pid, function(i) {
    f <- cfg$V1[cfg$V2 == "1" & as.integer(cfg$V3) == i]; b <- cfg$V1[cfg$V2 == "0" & as.integer(cfg$V3) == i]
    if (length(f) != 1L || length(b) != 1L || is.null(obs_cand)) return(NA_real_)
    hit <- which((obs_cand$species1 == f & obs_cand$species2 == b) | (obs_cand$species1 == b & obs_cand$species2 == f))
    if (length(hit)) obs_cand$pss_score[hit[1]] else NA_real_
  }, numeric(1))
  if (length(obs_pss) != target_pairs || anyNA(obs_pss)) {
    log_msg("WARN", "the PSS of the observed pairs could not be recovered; PSS matching is off")
  } else {
    pss_profile <- obs_pss
    log_msg("INFO", sprintf("PSS matching on: observed PSS (selection order) %s, tolerance %.0f%%",
                            paste(format(obs_pss, digits = 4), collapse = " "), 100 * match_pss_tol))
    # Consistency with the hypotheses table written by the observed run, when it sits next to the config.
    .hpf <- file.path(if (dir.exists(config.file)) config.file else dirname(config.file), "contrast_hypotheses_pairs.tsv")
    if (file.exists(.hpf)) {
      .hp1 <- read.delim(.hpf, stringsAsFactors = FALSE)
      .hp1 <- .hp1[.hp1$hypothesis_id == "H1", ]
      if (nrow(.hp1) == target_pairs && "pss_score" %in% names(.hp1) &&
          max(abs(.hp1$pss_score[order(.hp1$pair)] - obs_pss)) > 1e-6) {
        log_msg("WARN", "the recomputed PSS of the observed pairs differs from contrast_hypotheses_pairs.tsv (H1)")
      }
    }
  }
}

# Pool and capacity of the observed trait under its own selection rule, for the harvest audit.
obs_cap <- tryCatch(
  lean_draw_capacity(
    starting.values, D, pruned.tree, cov_bm, cov_ou, selected_model,
    if (use_ci) setNames(phenotype.df$ci_lb, phenotype.df$species)[names(starting.values)] else NULL,
    if (use_ci) setNames(phenotype.df$ci_ub, phenotype.df$species)[names(starting.values)] else NULL,
    pss_top_pct,
    if (trait_type == "ordinal") TRUE else if (trait_type == "continuous") FALSE else NULL,
    if (use_ci) setNames(phenotype.df$n_pop, phenotype.df$species)[names(starting.values)] else NULL),
  error = function(err) c(pool_gated = NA_real_, pool_ungated = NA_real_, capacity = NA_real_))

# ── Harvest ───────────────────────────────────────────────────────────────────
# Draws are permulations of the observed trait: each is graded by
# evaluate_lean_contrast_selection() into Tier 1, Tier 2 or a rejection. The pool is the
# first `number.of.cycles` Tier-1 draws, completed with Tier 2 only when Tier 1 falls short.
start.time <- Sys.time()
log_msg("START", sprintf("Harvesting pool (strategy: %s, pool_size: %d, max_tries: %d)",
                         selection.strategy, number.of.cycles, max_tries))

rec_ord    <- order(starting.values)
rec_value  <- unname(starting.values[rec_ord])
n_rec      <- length(rec_value)
rank_order <- seq_len(n_rec)
rec_lb <- rec_ub <- rec_n <- NULL
if (use_ci) {
  ci_map <- setNames(phenotype.df$ci_lb, phenotype.df$species)
  rec_lb <- unname(ci_map[names(starting.values)[rec_ord]])
  ci_map <- setNames(phenotype.df$ci_ub, phenotype.df$species)
  rec_ub <- unname(ci_map[names(starting.values)[rec_ord]])
  n_map <- setNames(phenotype.df$n_pop, phenotype.df$species)
  rec_n <- unname(n_map[names(starting.values)[rec_ord]])
}

tier1 <- vector("list", number.of.cycles)
tier2 <- vector("list", number.of.cycles)
n1 <- 0L; n2 <- 0L
total_draws <- 0L
reject_reasons <- character(0)
# Harvest audit: the worst PSS mismatch of every draw (NA where no complete set of pairs was found), the
# design-matching rounds, and the pool and capacity of a sample of the draws (written to permulation_harvest.tsv).
draw_mm  <- rep(NA_real_, 10000L)
dm_rows  <- list()
cap_rows <- list()
CAP_SAMPLE <- 200L
ordinal_arg <- if (trait_type == "ordinal") TRUE else if (trait_type == "continuous") FALSE else NULL

# Budget escalation is target-aware. `max_tries` is meant to abort a trait whose Dunn
# geometry makes a full pool essentially unreachable (low Tier-1 acceptance), not to cap a
# run that is merely draw-starved because the pool is large and the starting budget small.
# When the inner loop stalls short of the pool:
#   * estimate the Tier-1 acceptance rate from the draws so far;
#   * if it is healthy, project the budget needed to finish (+ headroom) and keep
#     going, bounded only by HARVEST_HARD_CAP x pool_size;
#   * if it is genuinely poor, bump modestly, count it, and give up after
#     MAX_LOWRATE_ESCALATIONS with a diagnostic that names the acceptance rate.
HARVEST_HARD_CAP        <- 50L    # x candidate target: absolute draw ceiling
MIN_VIABLE_TIER1_RATE   <- 0.05   # below this the trait cannot realistically fill the pool
MAX_LOWRATE_ESCALATIONS <- 3L
LOWRATE_FACTOR          <- 1.5
BUDGET_HEADROOM         <- 1.15   # over the point estimate of draws still needed

budget              <- max_tries
lowrate_escalations <- 0L
escalations         <- 0L         # kept for the "filled after N escalation(s)" log below

# ── FOP harvest helpers (design matching and FOP mirror) ──────────────────────
# `lean_fop_harvest` seeds its draws with the pipeline seed (restoring the caller's RNG
# state afterwards) and `evaluate_lean_contrast_selection` is RNG-free, so the harvest of
# a cycle is a pure function of its own inputs: counting it here and re-harvesting it in
# the FOP mirror below yields the same hypotheses.
fop_hypotheses <- function(e, label) {
  if (is.null(e$pvec) || is.null(e$fg) || is.null(e$bg)) return(NULL)
  tryCatch(
    lean_fop_harvest(
      trait_vec = e$pvec, D = D, target_pairs = target_pairs,
      tree = pruned.tree, cov_bm = cov_bm, cov_ou = cov_ou,
      selected_model = selected_model,
      ci_lb = e$ci_lb_draw, ci_ub = e$ci_ub_draw, n_vec = e$n_draw,
      top_pct = pss_top_pct, max_fop = max_fop, seed = seed_arg,
      ordinal = if (trait_type == "ordinal") TRUE
                else if (trait_type == "continuous") FALSE else NULL,
      canon_pairs = if (!is.null(e$canon)) e$canon
                    else data.frame(species1 = e$fg, species2 = e$bg, stringsAsFactors = FALSE)),
    error = function(err) { log_msg("WARN", sprintf("FOP harvest %s: %s", label, conditionMessage(err))); NULL })
}
# Hypotheses a cycle contributes; an empty harvest falls back to H1 only.
fop_count <- function(e) {
  hv <- fop_hypotheses(e, paste0("draw ", e$draw_id))
  if (is.null(hv)) 1L else max(1L, length(hv$hypotheses))
}

# ── Design matching ───────────────────────────────────────────────────────────
# The observed statistic is pooled over the observed FOP harvest (n_hyp_obs hypotheses).
# A null cycle with fewer hypotheses has fewer chances to re-detect a position (empirical
# p biased low); one with more has more chances (biased high). With match_fop, each
# candidate is harvested with the same max_fop (hence the same search budget) as the
# observed harvest; the pool keeps, in draw order, only the candidates whose harvest
# reaches n_hyp_obs, and the FOP mirror below keeps the top n_hyp_obs hypotheses of each
# kept cycle. When too few candidates qualify, the candidate target is extended according
# to the observed match rate and harvesting resumes. If the rate is below MIN_MATCH_RATE
# or the target reaches MATCH_HARD_CAP x pool size, the shortfall is filled with the
# remaining candidates of largest harvest and the gap is logged.
MIN_MATCH_RATE <- 0.05
MATCH_HARD_CAP <- 20L
pool_target <- number.of.cycles
nhyp_cache  <- integer(0)   # draw_id -> FOP hypothesis count

repeat {
  hard_cap <- max(as.numeric(max_tries), HARVEST_HARD_CAP * pool_target)
  repeat {
    while (n1 < pool_target && total_draws < budget) {
      total_draws <- total_draws + 1L

      sim_v <- simulatevec(starting.values, simulation_tree)
      sim_ord <- order(sim_v)
      
      pvec <- setNames(rec_value[order(sim_ord)], names(sim_v))

      ci_lb_draw <- NULL; ci_ub_draw <- NULL; n_draw <- NULL
      if (use_ci) {
        ci_lb_draw <- setNames(rec_lb[order(sim_ord)], names(sim_v))
        ci_ub_draw <- setNames(rec_ub[order(sim_ord)], names(sim_v))
        n_draw     <- setNames(rec_n[order(sim_ord)], names(sim_v))
      }

      e <- tryCatch(
        evaluate_lean_contrast_selection(
          trait_vec = pvec, D = D, target_pairs = target_pairs,
          tree = pruned.tree, cov_bm = cov_bm, cov_ou = cov_ou,
          selected_model = selected_model,
          ci_lb = ci_lb_draw, ci_ub = ci_ub_draw,
          top_pct = pss_top_pct, n_vec = n_draw,
          pss_profile = pss_profile, pss_tol = match_pss_tol,
          ordinal = if (trait_type == "ordinal") TRUE
                    else if (trait_type == "continuous") FALSE
                    else NULL
        ),
        error = function(err) list(tier = 0L, n_pairs = 0L, dunn_min = 0,
                                   n_below = NA_integer_, fg = NULL, bg = NULL,
                                   reason = paste0("error: ", conditionMessage(err)))
      )

      if (total_draws > length(draw_mm)) draw_mm <- c(draw_mm, rep(NA_real_, length(draw_mm)))
      if (!is.null(e$mismatch)) draw_mm[total_draws] <- e$mismatch
      if (total_draws <= CAP_SAMPLE) {
        cap_rows[[total_draws]] <- tryCatch(
          lean_draw_capacity(pvec, D, pruned.tree, cov_bm, cov_ou, selected_model, ci_lb_draw, ci_ub_draw,
                             pss_top_pct, ordinal_arg, n_draw),
          error = function(err) c(pool_gated = NA_real_, pool_ungated = NA_real_, capacity = NA_real_))
      }

      # Accepted cycles keep their permuted vector (and CI/n draws), so that the FOP mirror
      # can harvest alternative hypotheses around this exact labeling once the pool is assembled.
      if (fop_null && e$tier %in% c(1L, 2L)) {
        e$draw_id <- total_draws
        e$pvec <- pvec
        e$ci_lb_draw <- ci_lb_draw; e$ci_ub_draw <- ci_ub_draw; e$n_draw <- n_draw
      }

      if (e$tier == 1L) {
        n1 <- n1 + 1L; tier1[[n1]] <- e
      } else if (e$tier == 2L && n2 < pool_target) {
        n2 <- n2 + 1L; tier2[[n2]] <- e
      } else if (e$tier == 0L) {
        reject_reasons <- c(reject_reasons, e$reason)
      }

      if (total_draws %% 5000 == 0 || n1 == pool_target) {
        el <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
        log_msg("PROGRESS", sprintf("draws=%d/%d | Tier1=%d/%d | Tier2=%d | %.0f draws/s | %.1f min",
                                    total_draws, budget, n1, pool_target, n2,
                                    total_draws / el, el / 60))
      }
    }

    if (n1 >= pool_target) break

    tier1_rate <- (n1 + 1) / (total_draws + 1)

    if (tier1_rate >= MIN_VIABLE_TIER1_RATE) {
      # Draw-starved, not rejection-bound: the budget is extended toward the projected
      # finish. This does not count against the low-acceptance abort counter.
      proj       <- total_draws + ceiling((pool_target - n1) / tier1_rate * BUDGET_HEADROOM)
      new_budget <- min(max(proj, ceiling(budget * 1.25)), hard_cap)
      if (new_budget <= budget) break   # already at the hard cap and still short
      log_msg("INFO", sprintf(
        paste0("Pool %d/%d after %d draws at %.1f%% Tier-1 acceptance — draw-starved, not ",
               "rejection-bound; extending budget %d -> %d (hard cap %d)"),
        n1, pool_target, total_draws, 100 * tier1_rate,
        as.integer(budget), as.integer(new_budget), as.integer(hard_cap)))
      budget <- new_budget
    } else {
      if (lowrate_escalations >= MAX_LOWRATE_ESCALATIONS) break
      lowrate_escalations <- lowrate_escalations + 1L
      escalations         <- escalations + 1L
      new_budget <- min(ceiling(budget * LOWRATE_FACTOR), hard_cap)
      if (new_budget <= budget) break
      budget <- new_budget
      log_msg("WARN", sprintf(
        paste0("Pool not filled from Tier 1 (%d/%d) after %d draws at only %.1f%% Tier-1 ",
               "acceptance — escalating budget to %d (escalation %d of %d)"),
        n1, pool_target, total_draws, 100 * tier1_rate,
        as.integer(budget), lowrate_escalations, MAX_LOWRATE_ESCALATIONS))
    }
  }

  # ── Assemble the candidate pool: Tier 1 first, Tier 2 to fill a shortfall ───
  if (n1 >= pool_target) {
    cands <- tier1[seq_len(pool_target)]
    log_msg("INFO", sprintf("Pool filled entirely from Tier 1 (%d/%d) in %d draws%s",
                            pool_target, pool_target, total_draws,
                            if (escalations) sprintf(" after %d escalation(s)", escalations) else ""))
  } else if (n1 + n2 >= pool_target) {
    use2 <- pool_target - n1
    cands <- c(tier1[seq_len(n1)], tier2[seq_len(use2)])
    log_msg("WARN", sprintf(
      "Tier 1 exhausted after %d draws: topping up with Tier 2. %d Tier-1 + %d Tier-2 = %d/%d records",
      n1, use2, length(cands), pool_target))
  } else if (n1 + n2 >= number.of.cycles) {
    # A design-matching extension ran out of draws: keep every candidate drawn.
    cands <- c(tier1[seq_len(n1)], tier2[seq_len(n2)])
    log_msg("WARN", sprintf(
      "Design-matching extension stopped at %d of %d candidates after %d draws",
      length(cands), pool_target, total_draws))
  } else {
    tab <- sort(table(reject_reasons), decreasing = TRUE)
    top <- seq_len(min(3, length(tab)))
    final_rate <- (n1 + 1) / (total_draws + 1)
    # The message names the failure mode, so the remedy is clear from the log alone.
    diag <- if (final_rate >= MIN_VIABLE_TIER1_RATE) sprintf(
        paste0("Tier-1 acceptance was healthy (%.1f%%) — the run was DRAW-STARVED and hit ",
               "the hard cap (%d). Raise --max_tries / MAX_TRIES (>= ~%d for this pool) or ",
               "lower --caas_full_perms."),
        100 * final_rate, as.integer(hard_cap),
        as.integer(ceiling(number.of.cycles / final_rate * BUDGET_HEADROOM)))
      else sprintf(
        paste0("Tier-1 acceptance was only %.1f%% — this trait's Dunn geometry cannot ",
               "realistically fill a pool this size; lower --caas_full_perms or relax the ",
               "contrast-selection strategy."),
        100 * final_rate)
    stop(sprintf(
      paste0("Permulation pool could not be filled: %d Tier-1 + %d Tier-2 = %d of the requested %d ",
             "after %d draws and %d low-acceptance escalation(s) (final budget %d, hard cap %d).\n",
             "  %s\n",
             "  Top rejection reasons: %s"),
      n1, n2, n1 + n2, number.of.cycles, total_draws, escalations,
      as.integer(budget), as.integer(hard_cap), diag,
      if (length(tab)) paste(sprintf("%s (%d)", names(tab)[top], as.integer(tab)[top]), collapse = "; ") else "none recorded"))
  }

  if (!match_fop) {
    pool <- cands[seq_len(number.of.cycles)]
    break
  }

  # ── Design matching: count each new candidate's FOP harvest ─────────────────
  cand_ids <- vapply(cands, function(e) as.integer(e$draw_id), integer(1))
  new_i <- which(!(as.character(cand_ids) %in% names(nhyp_cache)))
  if (length(new_i)) {
    nh_new <- parallel::mclapply(cands[new_i], fop_count,
                                 mc.cores = max(1L, min(n_cpus, length(new_i))),
                                 mc.preschedule = TRUE)
    nh_new <- vapply(nh_new, function(x) if (is.numeric(x)) as.integer(x) else 1L, integer(1))
    nhyp_cache[as.character(cand_ids[new_i])] <- nh_new
  }
  nh_c <- unname(nhyp_cache[as.character(cand_ids)])
  qual <- which(nh_c >= n_hyp_obs)
  match_rate <- length(qual) / length(cands)
  log_msg("INFO", sprintf("Design matching: %d/%d candidates reach >= %d FOP hypotheses (%.1f%%)",
                          length(qual), length(cands), n_hyp_obs, 100 * match_rate))
  dm_rows[[length(dm_rows) + 1L]] <- c(target = pool_target, candidates = length(cands), reach = length(qual))

  if (length(qual) >= number.of.cycles) {
    pool <- cands[qual[seq_len(number.of.cycles)]]
    break
  }
  if (length(cands) < pool_target || match_rate < MIN_MATCH_RATE ||
      pool_target >= MATCH_HARD_CAP * number.of.cycles) {
    rest <- setdiff(seq_along(cands), qual)
    rest <- rest[order(-nh_c[rest], rest)]
    fill <- rest[seq_len(min(length(rest), number.of.cycles - length(qual)))]
    pool <- cands[sort(c(qual, fill))]
    log_msg("WARN", sprintf(
      paste0("Design matching incomplete: %d cycles reach the observed %d FOP hypotheses; ",
             "%d filled with the largest remaining harvests (min %d hypotheses); match rate %.1f%%"),
      length(qual), n_hyp_obs, length(fill),
      if (length(fill)) min(nh_c[fill]) else NA_integer_, 100 * match_rate))
    break
  }
  pool_target <- as.integer(min(
    ceiling(length(cands) + (number.of.cycles - length(qual)) / match_rate * BUDGET_HEADROOM),
    MATCH_HARD_CAP * number.of.cycles))
  log_msg("INFO", sprintf("Design matching: extending the candidate target to %d", pool_target))
}

# ── Write resample chunks + a tier/Dunn manifest ──────────────────────────────
# Linear row binding: do.call(rbind, <list of n 1-row data.frames>) is quadratic and
# becomes the wall-clock bottleneck for pools of ~1e5 rows (the FOP mirror below produces
# up to max_fop times as many). data.table::rbindlist is linear.
rbind_fast <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (!length(lst)) return(NULL)
  as.data.frame(data.table::rbindlist(lst, use.names = TRUE, fill = TRUE),
                stringsAsFactors = FALSE)
}

# ── Harvest summary ───────────────────────────────────────────────────────────
# permulation_harvest.tsv, long format (section, key, x, value): what the harvest tried and discarded, for the
# "Null harvest" tab of the scoring report. Aggregates only; the per-cycle design is in permulation_manifest.tsv.
.hv <- list()
.add <- function(section, key, x = NA_real_, value = NA_real_)
  .hv[[length(.hv) + 1L]] <<- data.frame(section = section, key = key, x = x, value = value, stringsAsFactors = FALSE)
.add("run", "draws", value = total_draws)
.add("run", "tier1", value = n1)
.add("run", "tier2", value = n2)
.add("run", "pool", value = length(pool))
.add("run", "escalations", value = escalations)
.add("settings", "match_pss", value = as.numeric(!is.null(pss_profile)))
.add("settings", "match_pss_tol", value = match_pss_tol)
.add("settings", "top_pct", value = pss_top_pct)
.add("settings", "n_hyp_obs", value = n_hyp_obs)
.add("settings", "target_pairs", value = target_pairs)
.rt <- table(reject_reasons)
for (i in seq_along(.rt)) .add("reject", names(.rt)[i], value = as.numeric(.rt[[i]]))
.mm <- draw_mm[seq_len(total_draws)]
if (!is.null(pss_profile)) {
  .add("tolerance", "complete", value = sum(!is.na(.mm)))
  for (t in sort(unique(c(seq(0.05, 1, by = 0.05), seq(1.25, 3, by = 0.25), match_pss_tol))))
    .add("tolerance", "complete_within_tol", x = t, value = sum(.mm <= log1p(t), na.rm = TRUE))
  for (i in seq_along(pss_profile)) .add("observed", "pss", x = i, value = pss_profile[i])
}
for (i in seq_along(dm_rows)) for (k in names(dm_rows[[i]])) .add("design", k, x = i, value = dm_rows[[i]][[k]])
for (i in seq_along(cap_rows)) if (!is.null(cap_rows[[i]])) for (k in names(cap_rows[[i]])) .add("capacity", k, x = i, value = cap_rows[[i]][[k]])
for (k in names(obs_cap)) .add("observed", k, value = obs_cap[[k]])
write.table(data.table::rbindlist(.hv), file.path(outdir, "permulation_harvest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Streamed writers: only the current chunk is held in memory, never the whole pool. Each
# resample file holds `chunk.size` cycles; the manifest is appended chunk by chunk to a
# single open connection, with its header written once.
file.counter <- 1L; chunk.start <- 1L
chunk_rows     <- list()
chunk_manifest <- list()

man_con <- file(file.path(outdir, "permulation_manifest.tsv"), "w")
writeLines(paste(c("cycle", "tier", "n_pairs", "dunn_min", "n_below", "mode",
                   "fg_values", "bg_values", "mean_distance", "mean_abs_diff", "mean_pss", "pss_mismatch"),
                 collapse = "\t"), man_con)

flush_chunk <- function() {
  fp <- file.path(outdir, sprintf("resample_%03d.tab", file.counter))
  write.table(rbind_fast(chunk_rows), file = fp, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(rbind_fast(chunk_manifest), file = man_con, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  file.counter   <<- file.counter + 1L
  chunk_rows      <<- list()
  chunk_manifest  <<- list()
}

fmt <- function(x) paste(format(x, digits = 15, trim = TRUE), collapse = ",")
for (b in seq_along(pool)) {
  e <- pool[[b]]
  chunk_rows[[length(chunk_rows) + 1L]] <- data.frame(
    cycle = paste0("b_", b),
    fg    = paste(e$fg, collapse = ","),
    bg    = paste(e$bg, collapse = ","),
    stringsAsFactors = FALSE
  )
  chunk_manifest[[length(chunk_manifest) + 1L]] <- data.frame(
    cycle = paste0("b_", b), tier = e$tier,
    n_pairs = e$n_pairs, dunn_min = e$dunn_min,
    n_below = e$n_below, mode = e$mode,
    fg_values = fmt(e$fg_values),
    bg_values = fmt(e$bg_values),
    mean_distance = e$mean_pd, mean_abs_diff = e$mean_df, mean_pss = e$mean_pss,
    pss_mismatch = if (is.null(e$mismatch)) NA_real_ else e$mismatch,
    stringsAsFactors = FALSE)
  if (b - chunk.start + 1L >= chunk.size || b == length(pool)) {
    flush_chunk()
    chunk.start <- b + 1L
  }
}
close(man_con)

# ── FOP mirror: per-cycle alternative-hypothesis harvest ──────────────────────
# Mirrors the observed FOP harvest (selection_algorithm.R::fop_pair_sel.f) for every
# accepted permulation cycle, so that the null is pooled over hypotheses like the observed
# data (fop_pool.py pools the hypotheses of a position). H1 is the already-accepted canonical
# contrast of the cycle; H2..Hn are Dunn-independent alternatives drawn from the same
# Voronoi domains, ranked by minimum and mean PSS and then by Dunn, and capped at max_fop.
#   fop_labelings.tab : "<cycle>~H<m>" \t fg_csv \t bg_csv   (labelings of the fanned discovery)
#   fop_pairs.tsv     : cycle, hypothesis_id, pair (domain), species1, species2, pss_score
#
# Parallel and streamed: the harvest of a cycle is a pure function of its own inputs (see
# the helpers above), so forking the loop changes only the interleaving of cycles, which is
# kept by consuming the worker results in pool order. The output therefore does not depend
# on the number of workers. Each batch of cycles is harvested with mclapply and its rows are
# appended to open connections, so peak memory is one batch, not the whole row set (up to
# max_fop x pool size).
if (fop_null) {
  FOP_BATCH <- 1000L
  n_workers <- max(1L, min(n_cpus, length(pool)))
  log_msg("START", sprintf(
    "FOP mirror harvest for %d accepted cycles (max_fop=%d, %d worker(s), batch=%d)",
    length(pool), max_fop, n_workers, FOP_BATCH))

  lab_path  <- file.path(outdir, "fop_labelings.tab")
  pair_path <- file.path(outdir, "fop_pairs.tsv")
  lab_con  <- file(lab_path, "w")
  pair_con <- file(pair_path, "w")
  writeLines(paste(c("cycle", "hypothesis_id", "pair", "species1", "species2",
                     "pss_score"), collapse = "\t"), pair_con)

  # One cycle -> its preformatted label rows, pair rows, and hypothesis count.
  fop_one <- function(b) {
    e <- pool[[b]]
    cyc <- paste0("b_", b)
    if (is.null(e$pvec) || is.null(e$fg) || is.null(e$bg))
      return(list(lab = NULL, pair = NULL, n_hyp = 0L))
    hv <- fop_hypotheses(e, cyc)
    # Design matching: a matched cycle carries exactly the observed number of hypotheses,
    # its top n_hyp_obs in harvest order (H1, then by minimum PSS).
    if (match_fop && !is.null(hv) && length(hv$hypotheses) > n_hyp_obs)
      hv$hypotheses <- hv$hypotheses[seq_len(n_hyp_obs)]
    if (is.null(hv) || length(hv$hypotheses) == 0L) {
      # Fallback to H1 only, so that the cycle still enters the fanned discovery
      return(list(
        lab = data.frame(cycle = paste0(cyc, "~H1"),
                         fg = paste(e$fg, collapse = ","),
                         bg = paste(e$bg, collapse = ","),
                         stringsAsFactors = FALSE),
        pair = NULL, n_hyp = 0L))
    }
    lab_l <- vector("list", length(hv$hypotheses))
    pair_l <- vector("list", length(hv$hypotheses))
    for (i in seq_along(hv$hypotheses)) {
      h_id <- names(hv$hypotheses)[i]
      hd <- hv$hypotheses[[h_id]]
      lab_l[[i]] <- data.frame(
        cycle = paste0(cyc, "~", h_id),
        fg = paste(hd$species1, collapse = ","),
        bg = paste(hd$species2, collapse = ","),
        stringsAsFactors = FALSE)
      pair_l[[i]] <- data.frame(
        cycle = cyc, hypothesis_id = h_id,
        pair = if ("cluster" %in% names(hd)) as.integer(hd$cluster) else seq_len(nrow(hd)),
        species1 = hd$species1, species2 = hd$species2,
        pss_score = if ("pss_score" %in% names(hd)) hd$pss_score else NA_real_,
        stringsAsFactors = FALSE)
    }
    list(lab = rbind_fast(lab_l), pair = rbind_fast(pair_l),
         n_hyp = length(hv$hypotheses))
  }

  # First cycle of each batch; seq(1, 0, by = ...) is an error, and an empty pool (N = 0) has no batch.
  batch_starts <- function(n, size) if (n > 0L) seq(1L, n, by = size) else integer(0)

  n_hyp_tot <- 0L
  any_lab   <- FALSE
  any_pair  <- FALSE
  for (start in batch_starts(length(pool), FOP_BATCH)) {
    block <- start:min(start + FOP_BATCH - 1L, length(pool))
    res <- parallel::mclapply(block, fop_one, mc.cores = min(n_workers, length(block)),
                              mc.preschedule = TRUE)
    lab_batch  <- vector("list", length(res))
    pair_batch <- vector("list", length(res))
    for (i in seq_along(res)) {
      r <- res[[i]]
      if (inherits(r, "try-error") || is.null(r)) {
        log_msg("WARN", sprintf("FOP worker for cycle b_%d failed: %s",
                                block[i], as.character(r)))
        next
      }
      lab_batch[[i]]  <- r$lab
      pair_batch[[i]] <- r$pair
      n_hyp_tot <- n_hyp_tot + r$n_hyp
    }
    lab_df  <- rbind_fast(lab_batch)
    pair_df <- rbind_fast(pair_batch)
    if (!is.null(lab_df)) {
      write.table(lab_df, file = lab_con, sep = "\t",
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
      any_lab <- TRUE
    }
    if (!is.null(pair_df)) {
      write.table(pair_df, file = pair_con, sep = "\t",
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
      any_pair <- TRUE
    }
    rm(res, lab_batch, pair_batch, lab_df, pair_df)
  }
  close(lab_con); close(pair_con)
  # The files exist only when they received rows.
  if (!any_lab)  unlink(lab_path)
  if (!any_pair) unlink(pair_path)

  log_msg("COMPLETE", sprintf("FOP mirror: %d cycles -> %d hypothesis labelings (mean %.1f/cycle) -> fop_labelings.tab",
                              length(pool), n_hyp_tot,
                              if (length(pool)) n_hyp_tot / length(pool) else 0))
}

# ── Summary ───────────────────────────────────────────────────────────────────
tiers <- vapply(pool, function(e) e$tier, integer(1))
dunns <- vapply(pool, function(e) e$dunn_min, numeric(1))
elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "mins"))
log_msg("COMPLETE", sprintf(
  "%d records (Tier1=%d, Tier2=%d) from %d draws | acceptance %.2f%% | median overall Dunn %.3f | %.2f min",
  length(pool), sum(tiers == 1L), sum(tiers == 2L), total_draws,
  100 * length(pool) / total_draws, median(dunns), elapsed))
